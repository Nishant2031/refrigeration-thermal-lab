"""Discretized moving-boundary refrigerant/sensible-fluid heat exchanger."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from CoolProp.CoolProp import PropsSI

from .heat_exchangers import lmtd


@dataclass(frozen=True)
class MovingBoundaryResult:
    mode: str
    heat_transfer_w: float
    required_ua_w_k: float
    secondary_outlet_c: float
    minimum_approach_k: float
    refrigerant_enthalpy_j_kg: np.ndarray
    refrigerant_temperature_c: np.ndarray
    secondary_temperature_c: np.ndarray
    cumulative_heat_w: np.ndarray


def moving_boundary_required_ua(
    *, mode: str, fluid: str, pressure_pa: float,
    refrigerant_h_in_j_kg: float, refrigerant_h_out_j_kg: float,
    refrigerant_mass_flow_kg_s: float, secondary_inlet_c: float,
    secondary_capacity_rate_w_k: float, segments: int = 80,
) -> MovingBoundaryResult:
    """Find UA required for a specified duty using counterflow segment LMTDs.

    CoolProp temperature at each pressure/enthalpy node naturally represents
    desuperheating, two-phase, subcooling, and superheating zones.
    """
    if mode not in {"evaporator", "condenser"}:
        raise ValueError("mode must be evaporator or condenser.")
    if segments < 10 or min(refrigerant_mass_flow_kg_s, secondary_capacity_rate_w_k) <= 0:
        raise ValueError("Use >=10 segments and positive flow/capacity rate.")
    hin, hout = refrigerant_h_in_j_kg, refrigerant_h_out_j_kg
    expected = hout > hin if mode == "evaporator" else hout < hin
    if not expected:
        raise ValueError("Refrigerant enthalpy direction does not match exchanger mode.")
    h = np.linspace(hin, hout, segments + 1)
    tref = np.array([PropsSI("T", "P", pressure_pa, "Hmass", float(x), fluid) - 273.15 for x in h])
    total_q = refrigerant_mass_flow_kg_s * abs(hout - hin)
    qcum = refrigerant_mass_flow_kg_s * np.abs(h - hin)
    if mode == "evaporator":
        secondary_out = secondary_inlet_c - total_q / secondary_capacity_rate_w_k
        # x=0 is refrigerant inlet and secondary outlet; counterflow.
        tsec = secondary_out + qcum / secondary_capacity_rate_w_k
        delta = tsec - tref
    else:
        secondary_out = secondary_inlet_c + total_q / secondary_capacity_rate_w_k
        # x=0 is refrigerant inlet and secondary outlet; counterflow.
        tsec = secondary_out - qcum / secondary_capacity_rate_w_k
        delta = tref - tsec
    if np.any(delta <= 0):
        raise ValueError(f"Temperature cross/pinch: minimum approach {delta.min():.4g} K.")
    dq = total_q / segments
    ua_segments = np.array([dq / lmtd(float(delta[i]), float(delta[i + 1])) for i in range(segments)])
    return MovingBoundaryResult(
        mode, total_q, float(ua_segments.sum()), secondary_out,
        float(delta.min()), h, tref, tsec, qcum,
    )

