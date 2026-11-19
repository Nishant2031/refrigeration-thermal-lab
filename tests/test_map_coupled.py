import pandas as pd

from thermal_model.data_workflow import fit_compressor_map
from thermal_model.map_coupled import MapCoupledInputs, solve_map_coupled_system


def test_map_closed_cycle_derives_positive_mass_flow_and_map_power():
    rows = []
    for evap in (-10, -5, 0, 5):
        for cond in (35, 40, 45, 50):
            for speed in (3000, 3600, 4200):
                rows.append({"evaporating_temperature_c": evap, "condensing_temperature_c": cond, "speed_rpm": speed,
                    "capacity_w": 5000 + 20*evap - 15*(cond-35) + .1*(speed-3600),
                    "power_w": 1100 - 5*evap + 14*(cond-35) + .02*(speed-3600)})
    fit = fit_compressor_map(pd.DataFrame(rows), .75, 3)
    r = solve_map_coupled_system(MapCoupledInputs(speed_rpm=3600, source_inlet_c=10, sink_inlet_c=28,
        evaporator_ua_w_k=950, condenser_ua_w_k=1200), fit)
    assert r.solver_success
    assert r.cycle.inputs.refrigerant_mass_flow_kg_s > 0
    assert r.map_electrical_power_w > 0
    assert r.map_cop > 0
