"""LMTD and effectiveness-NTU heat-exchanger models."""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log


def lmtd(delta_t_1_k: float, delta_t_2_k: float) -> float:
    """Log-mean temperature difference for two positive terminal differences."""
    if delta_t_1_k <= 0 or delta_t_2_k <= 0:
        raise ValueError("Terminal temperature differences must be positive (no temperature cross).")
    if abs(delta_t_1_k - delta_t_2_k) < 1e-10 * max(delta_t_1_k, delta_t_2_k):
        return 0.5 * (delta_t_1_k + delta_t_2_k)
    return (delta_t_1_k - delta_t_2_k) / log(delta_t_1_k / delta_t_2_k)


def effectiveness_ntu(ntu: float, capacity_ratio: float, arrangement: str = "counterflow") -> float:
    """Return sensible-sensible exchanger effectiveness.

    ``capacity_ratio`` is C_min/C_max. Supported arrangements are counterflow,
    parallel, and crossflow_both_unmixed.
    """
    if ntu < 0 or not 0 <= capacity_ratio <= 1:
        raise ValueError("Require NTU >= 0 and 0 <= C_r <= 1.")
    cr = capacity_ratio
    if ntu == 0:
        return 0.0
    if arrangement == "counterflow":
        if abs(1 - cr) < 1e-12:
            return ntu / (1 + ntu)
        a = exp(-ntu * (1 - cr))
        return (1 - a) / (1 - cr * a)
    if arrangement == "parallel":
        return (1 - exp(-ntu * (1 + cr))) / (1 + cr)
    if arrangement == "crossflow_both_unmixed":
        if cr < 1e-12:
            return 1 - exp(-ntu)
        # Common Kays-London correlation; continuous at small NTU.
        return 1 - exp((exp(-cr * ntu**0.78) - 1) * ntu**0.22 / cr)
    raise ValueError(f"Unknown arrangement: {arrangement}")


@dataclass(frozen=True)
class SensibleHXResult:
    heat_transfer_w: float
    hot_outlet_c: float
    cold_outlet_c: float
    effectiveness: float
    ntu: float
    capacity_ratio: float
    lmtd_k: float
    ua_from_lmtd_w_k: float
    energy_balance_error_w: float


def sensible_hx(
    hot_inlet_c: float, cold_inlet_c: float,
    hot_capacity_rate_w_k: float, cold_capacity_rate_w_k: float,
    ua_w_k: float, arrangement: str = "counterflow",
) -> SensibleHXResult:
    """Solve a two-stream sensible heat exchanger with the epsilon-NTU method."""
    if hot_inlet_c <= cold_inlet_c:
        raise ValueError("Hot inlet must be warmer than cold inlet.")
    if min(hot_capacity_rate_w_k, cold_capacity_rate_w_k, ua_w_k) <= 0:
        raise ValueError("Capacity rates and UA must be positive.")
    cmin = min(hot_capacity_rate_w_k, cold_capacity_rate_w_k)
    cmax = max(hot_capacity_rate_w_k, cold_capacity_rate_w_k)
    ntu = ua_w_k / cmin
    cr = cmin / cmax
    epsilon = effectiveness_ntu(ntu, cr, arrangement)
    qmax = cmin * (hot_inlet_c - cold_inlet_c)
    q = epsilon * qmax
    tho = hot_inlet_c - q / hot_capacity_rate_w_k
    tco = cold_inlet_c + q / cold_capacity_rate_w_k
    if arrangement == "counterflow":
        d1, d2 = hot_inlet_c - tco, tho - cold_inlet_c
    elif arrangement == "parallel":
        d1, d2 = hot_inlet_c - cold_inlet_c, tho - tco
    else:
        # A single terminal LMTD is not a valid crossflow design relation.
        d1 = d2 = float("nan")
    lmt = lmtd(d1, d2) if arrangement != "crossflow_both_unmixed" else float("nan")
    ua_check = q / lmt if arrangement != "crossflow_both_unmixed" else float("nan")
    return SensibleHXResult(q, tho, tco, epsilon, ntu, cr, lmt, ua_check, 0.0)


@dataclass(frozen=True)
class PhaseChangeResult:
    mode: str
    heat_transfer_w: float
    sensible_outlet_c: float
    effectiveness: float
    ntu: float
    lmtd_k: float
    ua_from_lmtd_w_k: float


def phase_change_hx(
    mode: str, saturation_temperature_c: float, sensible_inlet_c: float,
    sensible_capacity_rate_w_k: float, ua_w_k: float,
) -> PhaseChangeResult:
    """Phase-change refrigerant against one sensible stream.

    For an evaporator the sensible stream cools; for a condenser it heats.
    The refrigerant is assumed isothermal, so C_r -> 0 and epsilon=1-exp(-NTU).
    """
    if mode not in {"evaporator", "condenser"}:
        raise ValueError("mode must be 'evaporator' or 'condenser'.")
    if min(sensible_capacity_rate_w_k, ua_w_k) <= 0:
        raise ValueError("Capacity rate and UA must be positive.")
    driving = (sensible_inlet_c - saturation_temperature_c) if mode == "evaporator" else (saturation_temperature_c - sensible_inlet_c)
    if driving <= 0:
        raise ValueError("Inlet temperature is on the wrong side of saturation temperature.")
    ntu = ua_w_k / sensible_capacity_rate_w_k
    eps = 1 - exp(-ntu)
    q = eps * sensible_capacity_rate_w_k * driving
    sign = -1 if mode == "evaporator" else 1
    tout = sensible_inlet_c + sign * q / sensible_capacity_rate_w_k
    d_in = driving
    d_out = abs(tout - saturation_temperature_c)
    lmt = lmtd(d_in, d_out)
    return PhaseChangeResult(mode, q, tout, eps, ntu, lmt, q / lmt)

