import pandas as pd
import pytest

from thermal_model.data_workflow import fit_compressor_map, score_test_points, validate_compressor_map


def map_data():
    rows = []
    for e in (-10, -5, 0, 5):
        for c in (35, 40, 45, 50):
            for speed in (3000, 3600, 4200):
                rows.append({"evaporating_temperature_c": e, "condensing_temperature_c": c, "speed_rpm": speed, "capacity_w": 3500 + 70 * e - 35 * (c - 35) + .5 * speed, "power_w": 800 - 12 * e + 20 * (c - 35) + .07 * speed})
    return pd.DataFrame(rows)


def test_map_fit_and_heldout_score_known_surface():
    fit = fit_compressor_map(map_data(), .75, 19)
    capacity, power = fit.predict(-2, 42, 3600)
    assert capacity == pytest.approx(3500 + 70 * -2 - 35 * 7 + .5 * 3600, rel=1e-10)
    assert power == pytest.approx(800 - 12 * -2 + 20 * 7 + .07 * 3600, rel=1e-10)
    assert fit.capacity_validation.mape_percent < 1e-8


def test_map_schema_rejects_duplicates():
    bad = pd.concat([map_data().iloc[:14], map_data().iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate"):
        validate_compressor_map(bad)


def test_test_point_scores():
    table = pd.DataFrame({"test_id": ["A", "B"], "predicted_capacity_w": [1000, 2100], "measured_capacity_w": [1000, 2000], "predicted_power_w": [500, 630], "measured_power_w": [500, 600]})
    scores = score_test_points(table)
    assert scores["capacity"].mape_percent == pytest.approx(2.5)
    assert scores["power"].mape_percent == pytest.approx(2.5)
