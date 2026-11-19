"""Steady-state single-stage vapor-compression refrigeration cycle."""

from __future__ import annotations

from dataclasses import dataclass

from CoolProp.CoolProp import PropsSI

from .state import StatePoint


@dataclass(frozen=True)
class CycleInputs:
    fluid: str = "R134a"
    evaporating_temperature_c: float = -5.0
    condensing_temperature_c: float = 40.0
    superheat_k: float = 5.0
    subcooling_k: float = 3.0
    compressor_isentropic_efficiency: float = 0.72
    refrigerant_mass_flow_kg_s: float = 0.03

    def validate(self) -> None:
        if self.condensing_temperature_c <= self.evaporating_temperature_c:
            raise ValueError("Condensing temperature must exceed evaporating temperature.")
        if self.superheat_k < 0 or self.subcooling_k < 0:
            raise ValueError("Superheat and subcooling cannot be negative.")
        if self.superheat_k > 40 or self.subcooling_k > 40:
            raise ValueError("Superheat and subcooling above 40 K are outside this model's intended range.")
        if not 0 < self.compressor_isentropic_efficiency <= 1:
            raise ValueError("Compressor isentropic efficiency must be in (0, 1].")
        if self.refrigerant_mass_flow_kg_s <= 0:
            raise ValueError("Refrigerant mass flow must be positive.")


@dataclass(frozen=True)
class CycleResult:
    inputs: CycleInputs
    states: dict[str, StatePoint]
    state_2s: StatePoint
    evaporator_heat_w: float
    compressor_power_w: float
    condenser_heat_w: float
    cop_cooling: float
    cop_heating: float
    pressure_ratio: float
    energy_balance_error_w: float
    carnot_cop_cooling: float

    def summary(self) -> dict[str, float | str]:
        return {
            "fluid": self.inputs.fluid,
            "evaporating_temperature_C": self.inputs.evaporating_temperature_c,
            "condensing_temperature_C": self.inputs.condensing_temperature_c,
            "evaporator_heat_kW": self.evaporator_heat_w / 1000,
            "compressor_power_kW": self.compressor_power_w / 1000,
            "condenser_heat_kW": self.condenser_heat_w / 1000,
            "COP_cooling": self.cop_cooling,
            "COP_heating": self.cop_heating,
            "pressure_ratio": self.pressure_ratio,
            "energy_balance_error_W": self.energy_balance_error_w,
            "Carnot_COP_cooling": self.carnot_cop_cooling,
        }


def solve_cycle(inputs: CycleInputs) -> CycleResult:
    """Solve states 1-4 with constant-pressure heat exchangers and an isenthalpic valve.

    State numbering: 1 compressor inlet, 2 compressor outlet, 3 condenser
    outlet, 4 expansion-valve outlet. Pressure drops and heat loss are omitted.
    """
    inputs.validate()
    f = inputs.fluid
    te = inputs.evaporating_temperature_c + 273.15
    tc = inputs.condensing_temperature_c + 273.15
    # Saturation temperatures follow the compressor-rating convention
    # (AHRI 540 / EN 12900): both are DEW-point temperatures.  For pure fluids
    # dew and bubble coincide; for zeotropic blends such as R407C they differ
    # by the temperature glide, and subcooling is measured from the BUBBLE point.
    try:
        p_low = float(PropsSI("P", "T", te, "Q", 1, f))
        p_high = float(PropsSI("P", "T", tc, "Q", 1, f))
        t_bubble_high = float(PropsSI("T", "P", p_high, "Q", 0, f))
    except ValueError as exc:
        raise ValueError(
            f"{f} has no saturation state at {inputs.evaporating_temperature_c:.1f} °C / "
            f"{inputs.condensing_temperature_c:.1f} °C (check the critical temperature): {exc}") from exc
    if p_high <= p_low:
        raise ValueError("Computed high pressure is not above low pressure.")

    # A saturated state cannot be fixed by (P, T); use quality when superheat/subcooling is zero.
    if inputs.superheat_k > 1e-9:
        state1 = StatePoint.from_pair("1 compressor inlet", f, "P", p_low, "T", te + inputs.superheat_k)
    else:
        state1 = StatePoint.from_pair("1 compressor inlet", f, "P", p_low, "Q", 1.0)
    h2s = float(PropsSI("Hmass", "P", p_high, "Smass", state1.entropy_j_kg_k, f))
    h2 = state1.enthalpy_j_kg + (h2s - state1.enthalpy_j_kg) / inputs.compressor_isentropic_efficiency
    state2s = StatePoint.from_pair("2s isentropic outlet", f, "P", p_high, "Hmass", h2s)
    state2 = StatePoint.from_pair("2 compressor outlet", f, "P", p_high, "Hmass", h2)
    if inputs.subcooling_k > 1e-9:
        state3 = StatePoint.from_pair("3 condenser outlet", f, "P", p_high, "T", t_bubble_high - inputs.subcooling_k)
    else:
        state3 = StatePoint.from_pair("3 condenser outlet", f, "P", p_high, "Q", 0.0)
    state4 = StatePoint.from_pair("4 valve outlet", f, "P", p_low, "Hmass", state3.enthalpy_j_kg)

    mdot = inputs.refrigerant_mass_flow_kg_s
    q_e = mdot * (state1.enthalpy_j_kg - state4.enthalpy_j_kg)
    w_c = mdot * (state2.enthalpy_j_kg - state1.enthalpy_j_kg)
    q_c = mdot * (state2.enthalpy_j_kg - state3.enthalpy_j_kg)
    if min(q_e, w_c, q_c) <= 0:
        raise ValueError("Nonphysical cycle: heat duties and compressor work must be positive.")
    carnot = te / (tc - te)
    return CycleResult(
        inputs=inputs,
        states={"1": state1, "2": state2, "3": state3, "4": state4},
        state_2s=state2s,
        evaporator_heat_w=q_e,
        compressor_power_w=w_c,
        condenser_heat_w=q_c,
        cop_cooling=q_e / w_c,
        cop_heating=q_c / w_c,
        pressure_ratio=p_high / p_low,
        energy_balance_error_w=q_c - q_e - w_c,
        carnot_cop_cooling=carnot,
    )

