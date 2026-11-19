import math

import pytest

from thermal_model.heat_exchangers import effectiveness_ntu, lmtd, phase_change_hx, sensible_hx


def test_lmtd_known_value_and_equal_limit():
    assert lmtd(50.0, 30.0) == pytest.approx(39.1523, rel=1e-5)
    assert lmtd(20.0, 20.0) == pytest.approx(20.0)
    with pytest.raises(ValueError):
        lmtd(10.0, -1.0)


def test_counterflow_equal_capacity_limit():
    assert effectiveness_ntu(2.0, 1.0, "counterflow") == pytest.approx(2 / 3)


def test_sensible_hx_energy_and_lmtd_agree():
    r = sensible_hx(80, 20, 1800, 2400, 1100, "counterflow")
    qh = 1800 * (80 - r.hot_outlet_c)
    qc = 2400 * (r.cold_outlet_c - 20)
    assert qh == pytest.approx(qc, rel=1e-12)
    assert r.ua_from_lmtd_w_k == pytest.approx(1100, rel=1e-12)
    assert 0 < r.effectiveness < 1


@pytest.mark.parametrize("mode,tsat,tin", [("evaporator", -5, 12), ("condenser", 40, 30)])
def test_phase_change_epsilon_and_lmtd(mode, tsat, tin):
    r = phase_change_hx(mode, tsat, tin, 1000, 800)
    assert r.effectiveness == pytest.approx(1 - math.exp(-0.8))
    assert r.ua_from_lmtd_w_k == pytest.approx(800, rel=1e-12)
    assert abs(r.sensible_outlet_c - tsat) < abs(tin - tsat)

