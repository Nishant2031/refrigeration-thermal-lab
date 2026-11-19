"""Map-driven coupled cycle closure.

The manufacturer map closes refrigerant mass flow and electrical power.  The
finite-UA exchangers then determine the saturation temperatures.  Map use is
rejected outside the uploaded envelope - it is interpolation, not magic.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import inf

from .cycle import CycleInputs, CycleResult, solve_cycle
from .data_workflow import CompressorMapFit
from .moving_boundary import MovingBoundaryResult, moving_boundary_required_ua
from .solver import InfeasibleOperatingPoint, close_saturation_temperatures, cycle_temperature_bounds


@dataclass(frozen=True)
class MapCoupledInputs:
    fluid: str = "R134a"
    speed_rpm: float = 3600.0
    superheat_k: float = 5.0
    subcooling_k: float = 3.0
    source_inlet_c: float = 12.0
    source_capacity_rate_w_k: float = 1000.0
    evaporator_ua_w_k: float = 950.0
    sink_inlet_c: float = 30.0
    sink_capacity_rate_w_k: float = 1250.0
    condenser_ua_w_k: float = 1200.0
    auxiliary_power_w: float = 0.0
    segments: int = 80
    # Conditions at which the manufacturer map was rated.  Map capacity is
    # converted to mass flow with the RATING refrigerating effect (AHRI 540
    # practice); actual capacity then uses the OPERATING superheat/subcooling.
    # ``None`` means "same as the operating values".
    map_rating_superheat_k: float | None = None
    map_rating_subcooling_k: float | None = None


@dataclass(frozen=True)
class MapCoupledResult:
    inputs: MapCoupledInputs
    cycle: CycleResult
    evaporator: MovingBoundaryResult
    condenser: MovingBoundaryResult
    map_capacity_w: float
    map_electrical_power_w: float
    map_cop: float
    system_cop: float
    residual_ua_fraction: tuple[float, float]
    solver_success: bool
    message: str
    function_evaluations: int


def solve_map_coupled_system(inputs: MapCoupledInputs, compressor_map: CompressorMapFit) -> MapCoupledResult:
    """Solve a map-closed, fixed-speed steady operating point.

    The map's predicted evaporator capacity is converted to mass flow using the
    calculated refrigerating effect at each trial saturation condition.  Map
    electrical power is retained as the compressor electrical boundary.
    """
    if min(inputs.speed_rpm, inputs.source_capacity_rate_w_k, inputs.sink_capacity_rate_w_k,
           inputs.evaporator_ua_w_k, inputs.condenser_ua_w_k) <= 0:
        raise ValueError("Speed, UA values, and secondary capacity rates must be positive.")
    if inputs.auxiliary_power_w < 0:
        raise ValueError("Auxiliary power cannot be negative.")
    rating_sh = inputs.superheat_k if inputs.map_rating_superheat_k is None else inputs.map_rating_superheat_k
    rating_sc = inputs.subcooling_k if inputs.map_rating_subcooling_k is None else inputs.map_rating_subcooling_k
    if rating_sh < 0 or rating_sc < 0:
        raise ValueError("Map rating superheat and subcooling cannot be negative.")
    (te_lo, te_hi), (tc_lo, tc_hi) = cycle_temperature_bounds(inputs.fluid, inputs.source_inlet_c, inputs.sink_inlet_c)
    # Never extrapolate the map: intersect the physical window with its envelope.
    te_bounds = (max(te_lo, compressor_map.evaporating_range_c[0]), min(te_hi, compressor_map.evaporating_range_c[1]))
    tc_bounds = (max(tc_lo, compressor_map.condensing_range_c[0]), min(tc_hi, compressor_map.condensing_range_c[1]))
    if not compressor_map.speed_range_rpm[0] <= inputs.speed_rpm <= compressor_map.speed_range_rpm[1]:
        raise ValueError(f"Speed {inputs.speed_rpm:.0f} rpm is outside the map's {compressor_map.speed_range_rpm[0]:.0f}-"
                         f"{compressor_map.speed_range_rpm[1]:.0f} rpm range.")
    if not (te_bounds[0] < te_bounds[1] and tc_bounds[0] < tc_bounds[1]):
        raise InfeasibleOperatingPoint(
            "The compressor map envelope "
            f"(evaporating {compressor_map.evaporating_range_c[0]:.1f} to {compressor_map.evaporating_range_c[1]:.1f} °C, "
            f"condensing {compressor_map.condensing_range_c[0]:.1f} to {compressor_map.condensing_range_c[1]:.1f} °C) "
            "does not overlap the temperatures these source/sink conditions allow.")

    def evaluate(te: float, tc: float):
        q_map, w_map = compressor_map.predict(te, tc, inputs.speed_rpm)
        if q_map <= 0 or w_map <= 0:
            raise ValueError("Compressor map predicted nonphysical capacity or power.")
        rating = solve_cycle(CycleInputs(fluid=inputs.fluid, evaporating_temperature_c=te,
            condensing_temperature_c=tc, superheat_k=rating_sh,
            subcooling_k=rating_sc, refrigerant_mass_flow_kg_s=1.0))
        rating_effect = rating.states["1"].enthalpy_j_kg - rating.states["4"].enthalpy_j_kg
        if rating_effect <= 0:
            raise ValueError("Calculated refrigerating effect is nonphysical.")
        mdot = q_map / rating_effect
        cycle = solve_cycle(CycleInputs(fluid=inputs.fluid, evaporating_temperature_c=te,
            condensing_temperature_c=tc, superheat_k=inputs.superheat_k,
            subcooling_k=inputs.subcooling_k, refrigerant_mass_flow_kg_s=mdot))
        return cycle, mdot, q_map, w_map

    def exchanger(cycle, mdot, mode):
        if mode == "evaporator":
            return moving_boundary_required_ua(mode="evaporator", fluid=inputs.fluid,
                pressure_pa=cycle.states["1"].pressure_pa, refrigerant_h_in_j_kg=cycle.states["4"].enthalpy_j_kg,
                refrigerant_h_out_j_kg=cycle.states["1"].enthalpy_j_kg, refrigerant_mass_flow_kg_s=mdot,
                secondary_inlet_c=inputs.source_inlet_c, secondary_capacity_rate_w_k=inputs.source_capacity_rate_w_k,
                segments=inputs.segments)
        return moving_boundary_required_ua(mode="condenser", fluid=inputs.fluid,
            pressure_pa=cycle.states["2"].pressure_pa, refrigerant_h_in_j_kg=cycle.states["2"].enthalpy_j_kg,
            refrigerant_h_out_j_kg=cycle.states["3"].enthalpy_j_kg, refrigerant_mass_flow_kg_s=mdot,
            secondary_inlet_c=inputs.sink_inlet_c, secondary_capacity_rate_w_k=inputs.sink_capacity_rate_w_k,
            segments=inputs.segments)

    def required(te: float, tc: float) -> tuple[float, float]:
        try:
            cycle, mdot, _, _ = evaluate(te, tc)
        except (ValueError, OverflowError):
            return inf, inf
        out = []
        for mode in ("evaporator", "condenser"):
            try:
                out.append(exchanger(cycle, mdot, mode).required_ua_w_k)
            except (ValueError, OverflowError):
                out.append(inf)
        return out[0], out[1]

    closure = close_saturation_temperatures(
        required, evaporator_ua=inputs.evaporator_ua_w_k, condenser_ua=inputs.condenser_ua_w_k,
        te_bounds=te_bounds, tc_bounds=tc_bounds,
        initial=(inputs.source_inlet_c - 8.0, inputs.sink_inlet_c + 10.0))
    cycle, mdot, q_map, power = evaluate(closure.evaporating_c, closure.condensing_c)
    evap, cond = exchanger(cycle, mdot, "evaporator"), exchanger(cycle, mdot, "condenser")
    q = cycle.evaporator_heat_w  # capacity corrected to operating superheat/subcooling
    r = ((evap.required_ua_w_k-inputs.evaporator_ua_w_k)/inputs.evaporator_ua_w_k,
         (cond.required_ua_w_k-inputs.condenser_ua_w_k)/inputs.condenser_ua_w_k)
    success = max(map(abs, r)) < 1e-5
    return MapCoupledResult(inputs, cycle, evap, cond, q, power, q/power,
                            q/(power+inputs.auxiliary_power_w), r, success,
                            f"{closure.method}: {closure.message}", closure.function_evaluations)
