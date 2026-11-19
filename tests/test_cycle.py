import pytest

from thermal_model.cycle import CycleInputs, solve_cycle


def test_cycle_first_law_and_component_constraints():
    r = solve_cycle(CycleInputs())
    assert abs(r.energy_balance_error_w) < 1e-8 * r.condenser_heat_w
    assert r.states["4"].enthalpy_j_kg == pytest.approx(r.states["3"].enthalpy_j_kg, rel=1e-10)
    eta_recovered = ((r.state_2s.enthalpy_j_kg - r.states["1"].enthalpy_j_kg) /
                     (r.states["2"].enthalpy_j_kg - r.states["1"].enthalpy_j_kg))
    assert eta_recovered == pytest.approx(r.inputs.compressor_isentropic_efficiency, rel=1e-9)
    assert 0 < r.cop_cooling < r.carnot_cop_cooling
    assert r.cop_heating == pytest.approx(r.cop_cooling + 1, rel=1e-12)


def test_default_r134a_regression_window():
    """Broad physical range catches unit/state errors without pinning a CoolProp build."""
    r = solve_cycle(CycleInputs())
    assert 2.5 < r.cop_cooling < 5.5
    assert 2.5 < r.pressure_ratio < 5.5
    assert 4_000 < r.evaporator_heat_w < 7_000


@pytest.mark.parametrize("eta", [0, -0.2, 1.1])
def test_invalid_efficiency_rejected(eta):
    with pytest.raises(ValueError):
        solve_cycle(CycleInputs(compressor_isentropic_efficiency=eta))

