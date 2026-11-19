"""Run all model types and save numeric and graphical results."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from thermal_model import (
    CoupledInputs, CycleInputs, phase_change_hx, sensible_hx,
    solve_coupled_system, solve_cycle,
)
from thermal_model.plots import plot_cop_sweep, plot_cycle_ph, plot_hx_profiles


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    cycle = solve_cycle(CycleInputs())
    sensible = sensible_hx(80.0, 20.0, 1800.0, 2400.0, 1100.0, "counterflow")
    phase_evap = phase_change_hx("evaporator", -5.0, 12.0, 1000.0, 950.0)
    phase_cond = phase_change_hx("condenser", 40.0, 30.0, 1250.0, 1200.0)
    coupled = solve_coupled_system(CoupledInputs())
    if not coupled.solver_success:
        raise RuntimeError(f"Coupled solve failed: {coupled.solver_message}; residuals={coupled.residual_ua_fraction}")

    report = {
        "standalone_cycle": cycle.summary(),
        "sensible_hx": asdict(sensible),
        "phase_change_evaporator": asdict(phase_evap),
        "phase_change_condenser": asdict(phase_cond),
        "coupled_cycle": coupled.cycle.summary(),
        "coupled_heat_exchangers": {
            "evaporator_required_UA_W_K": coupled.evaporator.required_ua_w_k,
            "evaporator_available_UA_W_K": coupled.inputs.evaporator_ua_w_k,
            "source_outlet_C": coupled.evaporator.secondary_outlet_c,
            "evaporator_minimum_approach_K": coupled.evaporator.minimum_approach_k,
            "condenser_required_UA_W_K": coupled.condenser.required_ua_w_k,
            "condenser_available_UA_W_K": coupled.inputs.condenser_ua_w_k,
            "sink_outlet_C": coupled.condenser.secondary_outlet_c,
            "condenser_minimum_approach_K": coupled.condenser.minimum_approach_k,
            "solver_function_evaluations": coupled.function_evaluations,
            "fractional_UA_residuals": coupled.residual_ua_fraction,
        },
    }
    (RESULTS / "example_results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    with (RESULTS / "cycle_state_points.csv").open("w", newline="", encoding="utf-8") as f:
        rows = [cycle.states[k].as_display_dict() for k in ("1", "2", "3", "4")]
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    plot_cycle_ph(cycle, RESULTS / "cycle_ph_diagram.png")
    plot_hx_profiles(coupled, RESULTS / "coupled_hx_temperature_profiles.png")
    plot_cop_sweep(cycle.inputs, RESULTS / "cop_vs_condensing_temperature.png")

    print(json.dumps(report, indent=2))
    print(f"\nSaved outputs in {RESULTS}")


if __name__ == "__main__":
    main()

