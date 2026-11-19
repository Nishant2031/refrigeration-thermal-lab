"""Couple the refrigerant cycle to finite-UA source and sink heat exchangers."""

from __future__ import annotations

from dataclasses import dataclass

from math import inf

from .cycle import CycleInputs, CycleResult, solve_cycle
from .moving_boundary import MovingBoundaryResult, moving_boundary_required_ua
from .solver import InfeasibleOperatingPoint, close_saturation_temperatures, cycle_temperature_bounds


@dataclass(frozen=True)
class CoupledInputs:
    fluid: str = "R134a"
    refrigerant_mass_flow_kg_s: float = 0.03
    compressor_isentropic_efficiency: float = 0.72
    superheat_k: float = 5.0
    subcooling_k: float = 3.0
    source_inlet_c: float = 12.0
    source_capacity_rate_w_k: float = 1000.0
    evaporator_ua_w_k: float = 950.0
    sink_inlet_c: float = 30.0
    sink_capacity_rate_w_k: float = 1250.0
    condenser_ua_w_k: float = 1200.0
    auxiliary_power_w: float = 0.0
    case_name: str = "unnamed coupled steady-state case"
    model_mode: str = "known_ua_component_prediction"
    segments: int = 80


@dataclass(frozen=True)
class CoupledResult:
    inputs: CoupledInputs
    cycle: CycleResult
    evaporator: MovingBoundaryResult
    condenser: MovingBoundaryResult
    solver_success: bool
    solver_message: str
    residual_ua_fraction: tuple[float, float]
    function_evaluations: int

    @property
    def total_electrical_power_w(self) -> float:
        """Compressor plus explicitly declared fan/pump/control electrical power."""
        return self.cycle.compressor_power_w + self.inputs.auxiliary_power_w

    @property
    def system_cop_cooling(self) -> float:
        """Useful evaporator duty divided by the declared electrical boundary."""
        return self.cycle.evaporator_heat_w / self.total_electrical_power_w

    def traceable_summary(self) -> dict[str, object]:
        """A compact record suitable for an export or validation report."""
        return {
            "case_name": self.inputs.case_name,
            "model_mode": self.inputs.model_mode,
            "fluid": self.inputs.fluid,
            "solver_success": self.solver_success,
            "solver_message": self.solver_message,
            "function_evaluations": self.function_evaluations,
            "cooling_capacity_w": self.cycle.evaporator_heat_w,
            "compressor_power_w": self.cycle.compressor_power_w,
            "auxiliary_power_w": self.inputs.auxiliary_power_w,
            "total_electrical_power_w": self.total_electrical_power_w,
            "cycle_cop": self.cycle.cop_cooling,
            "system_cop": self.system_cop_cooling,
            "evaporator_minimum_approach_k": self.evaporator.minimum_approach_k,
            "condenser_minimum_approach_k": self.condenser.minimum_approach_k,
            "ua_residual_fraction": self.residual_ua_fraction,
        }


def _cycle(inp: CoupledInputs, te_c: float, tc_c: float) -> CycleResult:
    return solve_cycle(CycleInputs(
        fluid=inp.fluid,
        evaporating_temperature_c=te_c,
        condensing_temperature_c=tc_c,
        superheat_k=inp.superheat_k,
        subcooling_k=inp.subcooling_k,
        compressor_isentropic_efficiency=inp.compressor_isentropic_efficiency,
        refrigerant_mass_flow_kg_s=inp.refrigerant_mass_flow_kg_s,
    ))


def _exchanger(inp: CoupledInputs, cyc: CycleResult, mode: str) -> MovingBoundaryResult:
    if mode == "evaporator":
        return moving_boundary_required_ua(
            mode="evaporator", fluid=inp.fluid,
            pressure_pa=cyc.states["1"].pressure_pa,
            refrigerant_h_in_j_kg=cyc.states["4"].enthalpy_j_kg,
            refrigerant_h_out_j_kg=cyc.states["1"].enthalpy_j_kg,
            refrigerant_mass_flow_kg_s=inp.refrigerant_mass_flow_kg_s,
            secondary_inlet_c=inp.source_inlet_c,
            secondary_capacity_rate_w_k=inp.source_capacity_rate_w_k,
            segments=inp.segments,
        )
    return moving_boundary_required_ua(
        mode="condenser", fluid=inp.fluid,
        pressure_pa=cyc.states["2"].pressure_pa,
        refrigerant_h_in_j_kg=cyc.states["2"].enthalpy_j_kg,
        refrigerant_h_out_j_kg=cyc.states["3"].enthalpy_j_kg,
        refrigerant_mass_flow_kg_s=inp.refrigerant_mass_flow_kg_s,
        secondary_inlet_c=inp.sink_inlet_c,
        secondary_capacity_rate_w_k=inp.sink_capacity_rate_w_k,
        segments=inp.segments,
    )


def _exchangers(inp: CoupledInputs, cyc: CycleResult) -> tuple[MovingBoundaryResult, MovingBoundaryResult]:
    return _exchanger(inp, cyc, "evaporator"), _exchanger(inp, cyc, "condenser")


def _required_ua(inp: CoupledInputs, cyc: CycleResult) -> tuple[float, float]:
    """Required UA of each exchanger, evaluated independently; +inf marks a pinch or cross."""
    out = []
    for mode in ("evaporator", "condenser"):
        try:
            out.append(_exchanger(inp, cyc, mode).required_ua_w_k)
        except (ValueError, OverflowError):
            out.append(inf)
    return out[0], out[1]


def _evaluate(inp: CoupledInputs, te_c: float, tc_c: float):
    cyc = _cycle(inp, te_c, tc_c)
    ev, co = _exchangers(inp, cyc)
    return cyc, ev, co


def solve_coupled_system(inp: CoupledInputs = CoupledInputs()) -> CoupledResult:
    """Solve saturation temperatures so required and available UA match.

    Raises ``InfeasibleOperatingPoint`` (a ``ValueError``) with a specific
    explanation when no steady operating point exists.
    """
    if min(inp.evaporator_ua_w_k, inp.condenser_ua_w_k,
           inp.source_capacity_rate_w_k, inp.sink_capacity_rate_w_k,
           inp.refrigerant_mass_flow_kg_s) <= 0:
        raise ValueError("Refrigerant flow, UA values, and secondary capacity rates must be positive.")
    if inp.auxiliary_power_w < 0:
        raise ValueError("Auxiliary power cannot be negative.")
    if not 0 < inp.compressor_isentropic_efficiency <= 1:
        raise ValueError("Compressor isentropic efficiency must be in (0, 1].")
    te_bounds, tc_bounds = cycle_temperature_bounds(inp.fluid, inp.source_inlet_c, inp.sink_inlet_c)

    def required(te: float, tc: float) -> tuple[float, float]:
        try:
            cyc = _cycle(inp, te, tc)
        except (ValueError, OverflowError):
            return inf, inf
        return _required_ua(inp, cyc)

    closure = close_saturation_temperatures(
        required, evaporator_ua=inp.evaporator_ua_w_k, condenser_ua=inp.condenser_ua_w_k,
        te_bounds=te_bounds, tc_bounds=tc_bounds,
        initial=(inp.source_inlet_c - 8.0, inp.sink_inlet_c + 10.0))
    cyc, ev, co = _evaluate(inp, closure.evaporating_c, closure.condensing_c)
    r = ((ev.required_ua_w_k - inp.evaporator_ua_w_k) / inp.evaporator_ua_w_k,
         (co.required_ua_w_k - inp.condenser_ua_w_k) / inp.condenser_ua_w_k)
    success = max(abs(r[0]), abs(r[1])) < 1e-5
    message = f"{closure.method}: {closure.message}"
    return CoupledResult(inp, cyc, ev, co, success, message, r, closure.function_evaluations)
