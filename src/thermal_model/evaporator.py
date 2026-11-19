"""Bounded evaporator calculations with explicit physical assumptions.

These functions deliberately separate *required cooling* from *available
heat-transfer capacity*.  They are not a substitute for a geometry-specific,
two-phase pressure-drop/boiling model.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, log

from CoolProp.CoolProp import PropsSI


@dataclass(frozen=True)
class LiquidTaskResult:
    required_duty_w: float
    mass_flow_kg_s: float
    enthalpy_in_j_kg: float
    enthalpy_target_j_kg: float
    method: str


def liquid_cooling_task(*, fluid: str, mass_flow_kg_s: float, inlet_c: float,
                        target_outlet_c: float, pressure_pa: float = 101325.0,
                        cp_j_kg_k: float | None = None) -> LiquidTaskResult:
    """Calculate cooling *required* for a steady sensible liquid stream.

    A supplied ``cp`` is an explicit constant-property simplification.  With
    ``cp=None``, CoolProp liquid enthalpies are used at the declared pressure.
    """
    if mass_flow_kg_s <= 0:
        raise ValueError("Liquid mass flow must be positive.")
    if target_outlet_c >= inlet_c:
        raise ValueError("Cooling target must be below the liquid inlet temperature.")
    if cp_j_kg_k is not None:
        if cp_j_kg_k <= 0:
            raise ValueError("Specific heat must be positive.")
        q = mass_flow_kg_s * cp_j_kg_k * (inlet_c - target_outlet_c)
        return LiquidTaskResult(q, mass_flow_kg_s, cp_j_kg_k * inlet_c,
                                cp_j_kg_k * target_outlet_c, "constant_cp")
    try:
        h_in = float(PropsSI("Hmass", "T", inlet_c + 273.15, "P", pressure_pa, fluid))
        h_out = float(PropsSI("Hmass", "T", target_outlet_c + 273.15, "P", pressure_pa, fluid))
    except ValueError as exc:
        raise ValueError(f"CoolProp could not calculate {fluid} liquid properties: {exc}") from exc
    return LiquidTaskResult(mass_flow_kg_s * (h_in - h_out), mass_flow_kg_s,
                            h_in, h_out, "CoolProp enthalpy")


@dataclass(frozen=True)
class SaturationScreeningResult:
    ua_w_k: float
    capacity_rate_w_k: float
    saturation_temperature_c: float
    secondary_inlet_c: float
    secondary_outlet_c: float
    effectiveness: float
    duty_w: float
    lmtd_k: float
    ua_check_w_k: float


def constant_saturation_evaporator(*, ua_w_k: float, capacity_rate_w_k: float,
                                   secondary_inlet_c: float,
                                   saturation_temperature_c: float) -> SaturationScreeningResult:
    """UA screening for one sensible stream against an isothermal evaporator.

    Assumes a fixed refrigerant saturation temperature, no refrigerant pressure
    loss, and no independent prediction of refrigerant flow/superheat.
    """
    if ua_w_k <= 0 or capacity_rate_w_k <= 0:
        raise ValueError("UA and secondary capacity rate must be positive.")
    driving = secondary_inlet_c - saturation_temperature_c
    if driving <= 0:
        raise ValueError("Secondary inlet must be warmer than evaporation saturation temperature.")
    ntu = ua_w_k / capacity_rate_w_k
    eps = 1.0 - exp(-ntu)
    duty = eps * capacity_rate_w_k * driving
    outlet = secondary_inlet_c - duty / capacity_rate_w_k
    d_in, d_out = driving, outlet - saturation_temperature_c
    lmtd = d_in if abs(d_in - d_out) < 1e-12 else (d_in - d_out) / log(d_in / d_out)
    return SaturationScreeningResult(ua_w_k, capacity_rate_w_k,
        saturation_temperature_c, secondary_inlet_c, outlet, eps, duty, lmtd, duty / lmtd)


@dataclass(frozen=True)
class RefrigerantFixtureResult:
    saturation_pressure_pa: float
    inlet_enthalpy_j_kg: float
    outlet_enthalpy_j_kg: float
    duty_w: float
    superheat_k: float
    outlet_phase: str


def refrigerant_evaporator_states(*, fluid: str, saturation_temperature_c: float,
                                  inlet_quality: float, outlet_temperature_c: float,
                                  refrigerant_mass_flow_kg_s: float) -> RefrigerantFixtureResult:
    """Property-consistent evaporator inlet/outlet states for pure fluids.

    Assumes no refrigerant pressure drop, so the outlet pressure equals the
    saturation pressure.  At that pressure the outlet temperature cannot be
    below saturation (that would be subcooled liquid, i.e. a condenser), and
    an outlet exactly at saturation is reported as saturated vapor - a wet
    outlet cannot be identified from temperature alone.
    """
    if not 0 <= inlet_quality <= 1 or refrigerant_mass_flow_kg_s <= 0:
        raise ValueError("Quality must be 0..1 and refrigerant mass flow must be positive.")
    tsat_k = saturation_temperature_c + 273.15
    try:
        p = float(PropsSI("P", "T", tsat_k, "Q", 1, fluid))
        t_bubble = float(PropsSI("T", "P", p, "Q", 0, fluid))
    except ValueError as exc:
        raise ValueError(f"{fluid} has no saturation state at {saturation_temperature_c:.1f} °C: {exc}") from exc
    if abs(t_bubble - tsat_k) > 0.05:
        raise ValueError(f"{fluid} is a zeotropic blend (glide {tsat_k - t_bubble:.1f} K). This pure-fluid "
                         "helper does not support blends; use a dew/bubble-aware mixture model.")
    h_in = float(PropsSI("Hmass", "P", p, "Q", inlet_quality, fluid))
    superheat = outlet_temperature_c - saturation_temperature_c
    if superheat < -1e-6:
        raise ValueError(f"Outlet temperature {outlet_temperature_c:.2f} °C is below saturation "
                         f"({saturation_temperature_c:.2f} °C). Without pressure drop the evaporator outlet "
                         "must be at or above saturation.")
    if superheat < 0.01:
        h_out, phase, superheat = float(PropsSI("Hmass", "P", p, "Q", 1, fluid)), "saturated vapor", 0.0
    else:
        h_out, phase = float(PropsSI("Hmass", "P", p, "T", outlet_temperature_c + 273.15, fluid)), "superheated vapor"
    duty = refrigerant_mass_flow_kg_s * (h_out - h_in)
    if duty <= 0:
        raise ValueError("Non-physical state pair: outlet enthalpy is not above inlet enthalpy.")
    return RefrigerantFixtureResult(p, h_in, h_out, duty, superheat, phase)
