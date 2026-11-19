import pytest

from thermal_model.coupled import CoupledInputs, solve_coupled_system


def test_coupled_ua_closure_and_energy_balance():
    inp = CoupledInputs(segments=50)
    r = solve_coupled_system(inp)
    assert r.solver_success, r.solver_message
    assert r.evaporator.required_ua_w_k == pytest.approx(inp.evaporator_ua_w_k, rel=2e-5)
    assert r.condenser.required_ua_w_k == pytest.approx(inp.condenser_ua_w_k, rel=2e-5)
    assert r.evaporator.heat_transfer_w == pytest.approx(r.cycle.evaporator_heat_w, rel=1e-10)
    assert r.condenser.heat_transfer_w == pytest.approx(r.cycle.condenser_heat_w, rel=1e-10)
    assert r.evaporator.minimum_approach_k > 0
    assert r.condenser.minimum_approach_k > 0


def test_more_condenser_ua_reduces_condensing_temperature():
    low = solve_coupled_system(CoupledInputs(condenser_ua_w_k=900, segments=40))
    high = solve_coupled_system(CoupledInputs(condenser_ua_w_k=1600, segments=40))
    assert high.cycle.inputs.condensing_temperature_c < low.cycle.inputs.condensing_temperature_c


def test_system_cop_uses_declared_auxiliary_power_and_export_is_traceable():
    result = solve_coupled_system(CoupledInputs(auxiliary_power_w=180.0, case_name="audit-case", segments=40))
    assert result.system_cop_cooling < result.cycle.cop_cooling
    assert result.total_electrical_power_w == pytest.approx(result.cycle.compressor_power_w + 180.0)
    record = result.traceable_summary()
    assert record["case_name"] == "audit-case"
    assert record["auxiliary_power_w"] == 180.0

