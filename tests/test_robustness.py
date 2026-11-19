"""Regression tests for the robustness fixes and hand-checked reference values."""
from math import exp, log
from pathlib import Path

import pandas as pd
import pytest
from CoolProp.CoolProp import PropsSI

from thermal_model import (CoupledInputs, CycleInputs, constant_saturation_evaporator,
                           liquid_cooling_task, refrigerant_evaporator_states, solve_coupled_system, solve_cycle)
from thermal_model.advanced import single_phase_pressure_drop
from thermal_model.data_workflow import fit_compressor_map, validate_compressor_map
from thermal_model.map_coupled import MapCoupledInputs, solve_map_coupled_system
from thermal_model.plate_evaporator import (NU_LAMINAR_PARALLEL_PLATES, PlateEvaporatorGeometry, channel_nusselt,
                                            gnielinski_nusselt, plate_evaporator_screening)
from thermal_model.solver import InfeasibleOperatingPoint

EXAMPLE_MAP = Path(__file__).resolve().parents[1] / "examples" / "copeland_zr34k3e_tfd_map.csv"


# ---------------------------------------------------------------- coupled solver
@pytest.mark.parametrize("fluid", ["R134a", "R407C", "R290", "R32", "R1234yf"])
def test_every_offered_refrigerant_converges_on_default_inputs(fluid):
    r = solve_coupled_system(CoupledInputs(fluid=fluid))
    assert r.solver_success
    c = r.cycle
    assert abs(c.energy_balance_error_w) < 1e-6 * c.condenser_heat_w
    assert c.cop_cooling < c.carnot_cop_cooling
    # Energy conservation on the secondary streams.
    assert r.evaporator.secondary_outlet_c == pytest.approx(12.0 - c.evaporator_heat_w / 1000.0, abs=1e-9)
    assert r.condenser.secondary_outlet_c == pytest.approx(30.0 + c.condenser_heat_w / 1250.0, abs=1e-9)


def test_default_r134a_case_is_unchanged_by_the_refactor():
    r = solve_coupled_system(CoupledInputs())
    assert r.cycle.evaporator_heat_w == pytest.approx(4747, rel=2e-3)
    assert r.cycle.cop_cooling == pytest.approx(5.284, rel=2e-3)


def test_r290_previously_crashed_now_solves():
    r = solve_coupled_system(CoupledInputs(fluid="R290"))
    assert r.solver_success and r.cycle.evaporator_heat_w > 0


def test_undersized_evaporator_reports_a_clear_infeasibility():
    with pytest.raises(InfeasibleOperatingPoint, match="evaporator"):
        solve_coupled_system(CoupledInputs(evaporator_ua_w_k=60))


def test_undersized_condenser_reports_a_clear_infeasibility():
    with pytest.raises(InfeasibleOperatingPoint, match="condenser"):
        solve_coupled_system(CoupledInputs(condenser_ua_w_k=60, sink_capacity_rate_w_k=150))


def test_sink_above_critical_temperature_is_rejected_up_front():
    with pytest.raises(ValueError, match="critical"):
        solve_coupled_system(CoupledInputs(fluid="R32", sink_inlet_c=75))


# ---------------------------------------------------------------- cycle
def test_zero_superheat_and_subcooling_do_not_crash():
    r = solve_cycle(CycleInputs(superheat_k=0, subcooling_k=0))
    assert r.states["1"].quality == pytest.approx(1.0)
    assert r.states["3"].quality == pytest.approx(0.0)


def test_blend_uses_dew_point_pressures_and_bubble_point_subcooling():
    r = solve_cycle(CycleInputs(fluid="R407C", evaporating_temperature_c=0, condensing_temperature_c=40, subcooling_k=3))
    p_dew = PropsSI("P", "T", 313.15, "Q", 1, "R407C")
    assert r.states["2"].pressure_pa == pytest.approx(p_dew, rel=1e-9)
    t_bubble = PropsSI("T", "P", p_dew, "Q", 0, "R407C")
    assert r.states["3"].temperature_k == pytest.approx(t_bubble - 3, abs=1e-6)


# ---------------------------------------------------------------- evaporator panel
def test_ua_screening_matches_hand_calculation():
    s = constant_saturation_evaporator(ua_w_k=950, capacity_rate_w_k=1000, secondary_inlet_c=12, saturation_temperature_c=4)
    eps = 1 - exp(-0.95)
    assert s.effectiveness == pytest.approx(eps)
    assert s.duty_w == pytest.approx(eps * 1000 * 8)
    d_out = 12 - eps * 8 - 4
    assert s.lmtd_k == pytest.approx((8 - d_out) / log(8 / d_out))
    assert s.ua_check_w_k == pytest.approx(950)


def test_water_cooling_task_matches_m_cp_dt():
    t = liquid_cooling_task(fluid="Water", mass_flow_kg_s=1000 / 4180, inlet_c=12, target_outlet_c=7)
    assert t.required_duty_w == pytest.approx(1000 / 4180 * 4195 * 5, rel=5e-3)


def test_outlet_below_saturation_is_rejected_not_labelled_superheated():
    with pytest.raises(ValueError, match="below saturation"):
        refrigerant_evaporator_states(fluid="R134a", saturation_temperature_c=5, inlet_quality=0.2,
                                      outlet_temperature_c=0, refrigerant_mass_flow_kg_s=0.004)


def test_outlet_exactly_at_saturation_is_saturated_vapor():
    r = refrigerant_evaporator_states(fluid="R134a", saturation_temperature_c=5, inlet_quality=0.2,
                                      outlet_temperature_c=5, refrigerant_mass_flow_kg_s=0.004)
    assert r.outlet_phase == "saturated vapor" and r.superheat_k == 0 and r.duty_w > 0


def test_any_zeotropic_blend_is_rejected_by_the_pure_fluid_helper():
    with pytest.raises(ValueError, match="zeotropic"):
        refrigerant_evaporator_states(fluid="R404A", saturation_temperature_c=5, inlet_quality=0.2,
                                      outlet_temperature_c=10, refrigerant_mass_flow_kg_s=0.004)


# ---------------------------------------------------------------- plate correlations
def test_plate_channel_laminar_nusselt_is_parallel_plate_value_not_round_tube():
    nu, regime = channel_nusselt(500, 7)
    assert nu == pytest.approx(7.54) and regime == "laminar"


def test_gnielinski_reference_value():
    # Hand calculation: f = (0.790 ln 1e4 - 1.64)^-2, Pr = 7
    assert gnielinski_nusselt(1e4, 7) == pytest.approx(79.5, abs=0.1)


def test_channel_nusselt_is_continuous_across_regime_limits():
    assert channel_nusselt(2299.999, 7)[0] == pytest.approx(channel_nusselt(2300.001, 7)[0], abs=1e-3)
    assert channel_nusselt(9999.999, 7)[0] == pytest.approx(channel_nusselt(10000.001, 7)[0], abs=1e-3)


def test_laminar_friction_constant_for_parallel_plates():
    common = dict(mass_flow_kg_s=0.01, density_kg_m3=1000, viscosity_pa_s=1e-3,
                  hydraulic_diameter_m=0.005, flow_area_m2=3e-4, length_m=0.5)
    tube = single_phase_pressure_drop(**common)
    plates = single_phase_pressure_drop(**common, laminar_friction_constant=96)
    assert plates.friction_factor == pytest.approx(96 / tube.reynolds)
    assert plates.pressure_drop_pa == pytest.approx(1.5 * tube.pressure_drop_pa)


def test_default_plate_case_uses_corrected_h():
    r = plate_evaporator_screening(geometry=PlateEvaporatorGeometry(.5, .12, 30, .0025), secondary_fluid="Water",
                                   secondary_mass_flow_kg_s=.24, secondary_inlet_c=12, refrigerant_saturation_c=4,
                                   refrigerant_h_w_m2_k=1800)
    k = r.secondary_h_w_m2_k * 0.005 / NU_LAMINAR_PARALLEL_PLATES
    assert 0.55 < k < 0.60  # water thermal conductivity, W/m-K
    assert r.secondary_regime == "laminar"


def test_plate_glycol_secondary_fluid_runs():
    r = plate_evaporator_screening(geometry=PlateEvaporatorGeometry(.5, .12, 30, .0025), secondary_fluid="INCOMP::MEG-30%",
                                   secondary_mass_flow_kg_s=.3, secondary_inlet_c=12, refrigerant_saturation_c=2,
                                   refrigerant_h_w_m2_k=1800)
    assert r.screening.duty_w > 0


# ---------------------------------------------------------------- compressor map
def test_map_in_btu_per_hour_is_rejected_by_carnot_check():
    data = pd.read_csv(EXAMPLE_MAP)
    data["capacity_w"] = data["capacity_btu_h"]  # the original labelling error
    with pytest.raises(ValueError, match="Carnot"):
        validate_compressor_map(data)


def test_copeland_example_fits_and_closes_below_carnot():
    fit = fit_compressor_map(pd.read_csv(EXAMPLE_MAP), 0.75, 2026)
    assert fit.capacity_validation.mape_percent < 2
    assert fit.power_validation.mape_percent < 2
    r = solve_map_coupled_system(MapCoupledInputs(
        fluid="R407C", map_rating_superheat_k=11.1, map_rating_subcooling_k=8.3,
        source_inlet_c=12, source_capacity_rate_w_k=2500, evaporator_ua_w_k=2500,
        sink_inlet_c=24, sink_capacity_rate_w_k=3000, condenser_ua_w_k=3000), fit)
    assert r.solver_success
    assert 0.3 < r.map_cop / r.cycle.carnot_cop_cooling < 0.8
    lo, hi = fit.evaporating_range_c
    assert lo <= r.cycle.inputs.evaporating_temperature_c <= hi


def test_map_rating_equal_to_operating_returns_map_capacity_exactly():
    fit = fit_compressor_map(pd.read_csv(EXAMPLE_MAP), 0.75, 2026)
    r = solve_map_coupled_system(MapCoupledInputs(
        fluid="R407C", source_inlet_c=12, source_capacity_rate_w_k=2500, evaporator_ua_w_k=2500,
        sink_inlet_c=24, sink_capacity_rate_w_k=3000, condenser_ua_w_k=3000), fit)
    q_map, _ = fit.predict(r.cycle.inputs.evaporating_temperature_c, r.cycle.inputs.condensing_temperature_c, 3600)
    assert r.map_capacity_w == pytest.approx(q_map, rel=1e-9)


def test_map_envelope_that_cannot_overlap_is_reported():
    fit = fit_compressor_map(pd.read_csv(EXAMPLE_MAP), 0.75, 2026)
    with pytest.raises(InfeasibleOperatingPoint, match="envelope"):
        solve_map_coupled_system(MapCoupledInputs(fluid="R407C", sink_inlet_c=70), fit)


# ---------------------------------------------------------------- design sweep
def test_evaporator_ua_sweep_shows_monotone_cop_with_diminishing_returns():
    from thermal_model.sweep import diminishing_returns, sweep_coupled
    df = sweep_coupled(CoupledInputs(), "evaporator_ua_w_k", [300, 800, 1300, 1800, 2300, 2800])
    assert df["feasible"].all()
    assert df["cycle_cop"].is_monotonic_increasing
    assert df["evaporating_c"].is_monotonic_increasing
    assert "diminishing returns" in diminishing_returns(df, "evaporator_ua_w_k", "cycle_cop")


def test_sweep_keeps_infeasible_points_with_reason():
    from thermal_model.sweep import sweep_coupled
    df = sweep_coupled(CoupledInputs(), "evaporator_ua_w_k", [40, 950])
    assert list(df["feasible"]) == [False, True]
    assert "evaporator" in df.loc[0, "reason"]


def test_map_envelope_is_the_convex_hull_not_the_rectangle():
    fit = fit_compressor_map(pd.read_csv(EXAMPLE_MAP), 0.75, 2026)
    assert fit.in_envelope(5, 40, 3600)
    assert not fit.in_envelope(-23, 65, 3600)  # cold-evaporating / hot-condensing corner is not published
    with pytest.raises(ValueError, match="envelope"):
        fit.predict(-23, 65, 3600)


# ---------------------------------------------------------------- manufacturer datasheet
DATASHEET = EXAMPLE_MAP.parent / "copeland_zr34k3e_tfd_datasheet.csv"


def test_datasheet_transcription_is_internally_consistent():
    d = pd.read_csv(DATASHEET)
    assert len(d) == 57
    # EER printed to 0.1; capacity to 100 Btu/h and power to 10 W.
    assert ((d.capacity_btu_h / d.power_w - d.eer_btu_wh).abs() < 0.15).all()


def test_model_reproduces_manufacturer_mass_flow_and_isentropic_efficiency():
    from thermal_model.datasheet_validation import validate_against_datasheet
    v = validate_against_datasheet(pd.read_csv(DATASHEET), fluid="R407C",
                                   rating_superheat_k=20 / 1.8, rating_subcooling_k=15 / 1.8)
    assert v.mass_flow.mape_percent < 0.5
    assert v.mass_flow.maximum_absolute_percent_error < 1.0
    assert v.table.isentropic_eff_error_points.abs().max() < 0.5


def test_wrong_rating_conditions_are_detectable_in_the_datasheet_check():
    from thermal_model.datasheet_validation import validate_against_datasheet
    v = validate_against_datasheet(pd.read_csv(DATASHEET), fluid="R407C", rating_superheat_k=5, rating_subcooling_k=3)
    assert v.mass_flow.mape_percent > 5
