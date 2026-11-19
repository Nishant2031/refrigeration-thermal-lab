"""Traceable compressor-map fitting and test-point validation utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


COMPRESSOR_MAP_COLUMNS = ("evaporating_temperature_c", "condensing_temperature_c", "speed_rpm", "capacity_w", "power_w")
TEST_POINT_COLUMNS = ("test_id", "predicted_capacity_w", "measured_capacity_w", "predicted_power_w", "measured_power_w")


def _numeric(data: pd.DataFrame, columns: tuple[str, ...], label: str) -> pd.DataFrame:
    missing = [c for c in columns if c not in data.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
    frame = data.copy()
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if len(frame) == 0 or frame[list(columns)].isna().any().any():
        raise ValueError(f"{label} has no records or contains blank/nonnumeric required values.")
    return frame


def validate_compressor_map(data: pd.DataFrame) -> pd.DataFrame:
    frame = _numeric(data, COMPRESSOR_MAP_COLUMNS, "Compressor map")
    if (frame[["capacity_w", "power_w", "speed_rpm"]] <= 0).any().any():
        raise ValueError("Capacity, power, and speed must be positive.")
    if (frame["condensing_temperature_c"] <= frame["evaporating_temperature_c"]).any():
        raise ValueError("Each map row needs condensing temperature above evaporating temperature.")
    if frame.duplicated(list(COMPRESSOR_MAP_COLUMNS[:3])).any():
        raise ValueError("Duplicate map points have identical evaporating temperature, condensing temperature, and speed.")
    # Second-law check: no real compressor can beat the Carnot cooling COP.
    # Violations almost always mean capacity in Btu/h (or kW) labelled as W.
    te_k = frame["evaporating_temperature_c"] + 273.15
    tc_k = frame["condensing_temperature_c"] + 273.15
    cop = frame["capacity_w"] / frame["power_w"]
    carnot = te_k / (tc_k - te_k)
    bad = frame.index[cop > carnot]
    if len(bad):
        i = bad[0]
        raise ValueError(
            f"Row {int(i) + 2}: capacity/power gives COP {cop[i]:.1f}, above the Carnot limit {carnot[i]:.1f} "
            f"at {frame.at[i, 'evaporating_temperature_c']:.1f} / {frame.at[i, 'condensing_temperature_c']:.1f} °C. "
            "Check units: capacity must be in W (1 Btu/h = 0.2931 W) and power in W.")
    return frame


def _features(frame: pd.DataFrame) -> np.ndarray:
    e = frame["evaporating_temperature_c"].to_numpy(float)
    c = frame["condensing_temperature_c"].to_numpy(float)
    n = frame["speed_rpm"].to_numpy(float) / 1000.0
    return np.column_stack((np.ones(len(frame)), e, c, n, e * e, c * c, n * n, e * c, e * n, c * n))


@dataclass(frozen=True)
class ErrorMetrics:
    mae: float
    rmse: float
    mape_percent: float
    bias_percent: float
    maximum_absolute_percent_error: float


def error_metrics(predicted: np.ndarray, measured: np.ndarray) -> ErrorMetrics:
    predicted, measured = np.asarray(predicted, float), np.asarray(measured, float)
    if len(predicted) != len(measured) or len(measured) == 0 or np.any(measured == 0):
        raise ValueError("Predicted/measured values must have equal nonzero length; measured values cannot be zero.")
    error = predicted - measured
    percent = 100.0 * error / measured
    return ErrorMetrics(float(np.mean(np.abs(error))), float(np.sqrt(np.mean(error**2))), float(np.mean(np.abs(percent))), float(np.mean(percent)), float(np.max(np.abs(percent))))


@dataclass(frozen=True)
class CompressorMapFit:
    capacity_coefficients: np.ndarray
    power_coefficients: np.ndarray
    calibration_rows: int
    validation_rows: int
    capacity_validation: ErrorMetrics
    power_validation: ErrorMetrics
    random_seed: int
    evaporating_range_c: tuple[float, float]
    condensing_range_c: tuple[float, float]
    speed_range_rpm: tuple[float, float]
    # (evaporating, condensing) points of the map; used for a convex-hull envelope
    # check because real map envelopes are not rectangular.
    envelope_points: np.ndarray | None = None

    def in_envelope(self, evaporating_temperature_c: float, condensing_temperature_c: float, speed_rpm: float) -> bool:
        if not (self.evaporating_range_c[0] <= evaporating_temperature_c <= self.evaporating_range_c[1]
                and self.condensing_range_c[0] <= condensing_temperature_c <= self.condensing_range_c[1]
                and self.speed_range_rpm[0] <= speed_rpm <= self.speed_range_rpm[1]):
            return False
        if self.envelope_points is None or len(self.envelope_points) < 3:
            return True
        from scipy.spatial import Delaunay
        try:
            hull = Delaunay(self.envelope_points)
        except Exception:  # degenerate (collinear) map - fall back to the rectangle
            return True
        return bool(hull.find_simplex([[evaporating_temperature_c, condensing_temperature_c]], tol=1e-9)[0] >= 0)

    def predict(self, evaporating_temperature_c: float, condensing_temperature_c: float, speed_rpm: float) -> tuple[float, float]:
        if not self.in_envelope(evaporating_temperature_c, condensing_temperature_c, speed_rpm):
            raise ValueError("Requested point is outside the imported compressor-map envelope; extrapolation is disabled.")
        row = pd.DataFrame([{"evaporating_temperature_c": evaporating_temperature_c, "condensing_temperature_c": condensing_temperature_c, "speed_rpm": speed_rpm}])
        features = _features(row)[0]
        return float(features @ self.capacity_coefficients), float(features @ self.power_coefficients)


def fit_compressor_map(data: pd.DataFrame, calibration_fraction: float = 0.75, random_seed: int = 2026) -> CompressorMapFit:
    """Fit a transparent second-order map; score only held-out map rows."""
    frame = validate_compressor_map(data).reset_index(drop=True)
    if not 0.5 <= calibration_fraction < 1.0:
        raise ValueError("Calibration fraction must be in [0.5, 1).")
    if len(frame) < 14:
        raise ValueError("At least 14 distinct map rows are required for this second-order fit plus hold-out validation.")
    order = np.random.default_rng(random_seed).permutation(len(frame))
    split = min(max(10, int(round(len(frame) * calibration_fraction))), len(frame) - 1)
    cal, val = frame.iloc[order[:split]], frame.iloc[order[split:]]
    x_cal, x_val = _features(cal), _features(val)
    cap, *_ = np.linalg.lstsq(x_cal, cal["capacity_w"].to_numpy(float), rcond=None)
    power, *_ = np.linalg.lstsq(x_cal, cal["power_w"].to_numpy(float), rcond=None)
    return CompressorMapFit(cap, power, len(cal), len(val), error_metrics(x_val @ cap, val["capacity_w"].to_numpy(float)), error_metrics(x_val @ power, val["power_w"].to_numpy(float)), random_seed,
        (float(frame.evaporating_temperature_c.min()), float(frame.evaporating_temperature_c.max())),
        (float(frame.condensing_temperature_c.min()), float(frame.condensing_temperature_c.max())),
        (float(frame.speed_rpm.min()), float(frame.speed_rpm.max())),
        frame[["evaporating_temperature_c", "condensing_temperature_c"]].drop_duplicates().to_numpy(float))


def score_test_points(data: pd.DataFrame) -> dict[str, ErrorMetrics]:
    missing = [c for c in TEST_POINT_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(f"Test comparison is missing required columns: {', '.join(missing)}")
    if data["test_id"].astype(str).str.strip().eq("").any() or data["test_id"].duplicated().any():
        raise ValueError("test_id must be present and unique.")
    frame = _numeric(data, TEST_POINT_COLUMNS[1:], "Test comparison")
    if (frame[list(TEST_POINT_COLUMNS[1:])] <= 0).any().any():
        raise ValueError("Predicted and measured capacity/power values must be positive.")
    return {"capacity": error_metrics(frame["predicted_capacity_w"], frame["measured_capacity_w"]), "power": error_metrics(frame["predicted_power_w"], frame["measured_power_w"])}

