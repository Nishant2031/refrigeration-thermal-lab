import pytest

from thermal_model.evaporator import (constant_saturation_evaporator,
    liquid_cooling_task, refrigerant_evaporator_states)


def test_liquid_task_specification_fixture():
    mdot = 1000 / 4180
    assert liquid_cooling_task(fluid="Water", mass_flow_kg_s=mdot, inlet_c=12, target_outlet_c=10, cp_j_kg_k=4180).required_duty_w == pytest.approx(2000)
    assert liquid_cooling_task(fluid="Water", mass_flow_kg_s=mdot, inlet_c=12, target_outlet_c=7, cp_j_kg_k=4180).required_duty_w == pytest.approx(5000)
    assert liquid_cooling_task(fluid="Water", mass_flow_kg_s=mdot, inlet_c=12, target_outlet_c=5, cp_j_kg_k=4180).required_duty_w == pytest.approx(7000)


def test_constant_saturation_ua_fixture():
    r = constant_saturation_evaporator(ua_w_k=950, capacity_rate_w_k=1000, secondary_inlet_c=12, saturation_temperature_c=4)
    assert r.effectiveness == pytest.approx(0.6132589765, rel=1e-9)
    assert r.duty_w == pytest.approx(4906.071812, rel=1e-8)
    assert r.secondary_outlet_c == pytest.approx(7.0939281876, rel=1e-8)
    assert r.ua_check_w_k == pytest.approx(950, rel=1e-12)


def test_r134a_property_fixture():
    r = refrigerant_evaporator_states(fluid="R134a", saturation_temperature_c=5, inlet_quality=.2, outlet_temperature_c=10, refrigerant_mass_flow_kg_s=.004)
    assert r.saturation_pressure_pa == pytest.approx(349658.607861, rel=1e-6)
    assert r.inlet_enthalpy_j_kg == pytest.approx(245700.171530, rel=1e-6)
    assert r.outlet_enthalpy_j_kg == pytest.approx(406070.703967, rel=1e-6)
    assert r.duty_w == pytest.approx(641.482129749, rel=1e-6)
    assert r.superheat_k == pytest.approx(5)


def test_invalid_task_is_not_silently_zero():
    with pytest.raises(ValueError):
        liquid_cooling_task(fluid="Water", mass_flow_kg_s=0, inlet_c=12, target_outlet_c=7)
