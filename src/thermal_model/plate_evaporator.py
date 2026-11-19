"""Geometry-informed plate-evaporator screening path.

This is deliberately a *single supported preliminary geometry*: a counterflow
plate channel model with a prescribed refrigerant-side boiling coefficient.
It calculates secondary-side convection and pressure drop, resistance-network
UA, and a constant-saturation duty.  It does not claim a validated boiling or
two-phase pressure-drop correlation.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import log
from CoolProp.CoolProp import PropsSI
from .advanced import single_phase_pressure_drop
from .evaporator import constant_saturation_evaporator, SaturationScreeningResult


NU_LAMINAR_PARALLEL_PLATES = 7.54  # both walls at uniform temperature (Shah & London)
RE_LAMINAR_LIMIT = 2300.0
RE_TURBULENT_LIMIT = 1.0e4


def gnielinski_nusselt(reynolds: float, prandtl: float) -> float:
    """Gnielinski (1976) turbulent Nu with the Petukhov smooth-wall friction factor.

    Valid for 3000 <= Re <= 5e6 and 0.5 <= Pr <= 2000.
    """
    f = (0.790 * log(reynolds) - 1.64) ** -2
    return (f / 8) * (reynolds - 1000) * prandtl / (1 + 12.7 * (f / 8) ** 0.5 * (prandtl ** (2 / 3) - 1))


def channel_nusselt(reynolds: float, prandtl: float) -> tuple[float, str]:
    """Fully developed smooth parallel-plate channel Nusselt number and regime.

    Laminar: Nu = 7.54 (wide parallel plates, both walls at uniform
    temperature - the round-tube value 4.36 does NOT apply to a plate channel).
    Turbulent: Gnielinski.  Transition: linear interpolation in Re between the
    laminar value at Re = 2300 and Gnielinski at Re = 1e4 (Gnielinski, 2013).
    """
    if reynolds <= 0 or prandtl <= 0:
        raise ValueError("Reynolds and Prandtl numbers must be positive.")
    if reynolds < RE_LAMINAR_LIMIT:
        return NU_LAMINAR_PARALLEL_PLATES, "laminar"
    if reynolds >= RE_TURBULENT_LIMIT:
        return gnielinski_nusselt(reynolds, prandtl), "turbulent"
    gamma = (reynolds - RE_LAMINAR_LIMIT) / (RE_TURBULENT_LIMIT - RE_LAMINAR_LIMIT)
    nu = (1 - gamma) * NU_LAMINAR_PARALLEL_PLATES + gamma * gnielinski_nusselt(RE_TURBULENT_LIMIT, prandtl)
    return nu, "transitional"


@dataclass(frozen=True)
class PlateEvaporatorGeometry:
    plate_length_m: float
    plate_width_m: float
    plate_count: int
    channel_gap_m: float
    enlargement_factor: float = 1.15
    wall_thickness_m: float = 0.0006
    wall_conductivity_w_m_k: float = 15.0  # stainless steel
    passes: int = 1

    @property
    def channel_count(self) -> int:
        return self.plate_count - 1

    @property
    def secondary_channels(self) -> int:
        return self.channel_count // 2

    @property
    def active_area_m2(self) -> float:
        return (self.plate_count - 2) * self.plate_length_m * self.plate_width_m * self.enlargement_factor


@dataclass(frozen=True)
class PlateEvaporatorResult:
    secondary_h_w_m2_k: float
    secondary_reynolds: float
    secondary_nusselt: float
    secondary_pressure_drop_pa: float
    secondary_regime: str
    ua_w_k: float
    screening: SaturationScreeningResult
    correlation_status: str


def plate_evaporator_screening(*, geometry: PlateEvaporatorGeometry, secondary_fluid: str,
                               secondary_mass_flow_kg_s: float, secondary_inlet_c: float,
                               refrigerant_saturation_c: float, refrigerant_h_w_m2_k: float,
                               fouling_resistance_m2k_w: float = 0.0,
                               pressure_pa: float = 101325.0) -> PlateEvaporatorResult:
    if min(geometry.plate_length_m, geometry.plate_width_m, geometry.channel_gap_m,
           geometry.enlargement_factor, geometry.wall_thickness_m,
           geometry.wall_conductivity_w_m_k, geometry.passes,
           secondary_mass_flow_kg_s, refrigerant_h_w_m2_k) <= 0:
        raise ValueError("Geometry, mass flow, and refrigerant coefficient must be positive.")
    if int(geometry.plate_count) != geometry.plate_count:
        raise ValueError("Plate count must be a whole number.")
    if geometry.channel_gap_m > 0.02 or geometry.enlargement_factor > 2.0:
        raise ValueError("Channel gap above 20 mm or area enlargement above 2.0 is outside the plate-exchanger range of this model.")
    if geometry.plate_count < 4 or geometry.secondary_channels < 1:
        raise ValueError("At least four plates are needed to allocate refrigerant and secondary channels.")
    if fouling_resistance_m2k_w < 0:
        raise ValueError("Fouling resistance cannot be negative.")
    if secondary_inlet_c <= refrigerant_saturation_c:
        raise ValueError("Secondary inlet must be warmer than the refrigerant saturation temperature.")
    # Approximate bulk-mean temperature for properties (inlet and an outlet no colder than inlet - 10 K).
    t_k = 0.5 * (secondary_inlet_c + max(refrigerant_saturation_c, secondary_inlet_c - 10.0)) + 273.15
    try:
        rho = float(PropsSI("Dmass", "T", t_k, "P", pressure_pa, secondary_fluid))
        mu = float(PropsSI("VISCOSITY", "T", t_k, "P", pressure_pa, secondary_fluid))
        k = float(PropsSI("CONDUCTIVITY", "T", t_k, "P", pressure_pa, secondary_fluid))
        cp = float(PropsSI("Cpmass", "T", t_k, "P", pressure_pa, secondary_fluid))
    except ValueError as exc:
        raise ValueError(f"CoolProp could not evaluate {secondary_fluid} at {t_k - 273.15:.1f} °C "
                         "(it may be frozen or outside its valid range).") from exc
    # Wide parallel rectangular channels: Dh ≈ 2*gap and one channel area W*gap.
    dh, flow_area = 2 * geometry.channel_gap_m, geometry.plate_width_m * geometry.channel_gap_m
    mdot_channel = secondary_mass_flow_kg_s / geometry.secondary_channels
    dp = single_phase_pressure_drop(mass_flow_kg_s=mdot_channel, density_kg_m3=rho,
        viscosity_pa_s=mu, hydraulic_diameter_m=dh, flow_area_m2=flow_area,
        length_m=geometry.plate_length_m * geometry.passes, roughness_m=1.5e-6,
        laminar_friction_constant=96.0)  # parallel plates, not 64 (round tube)
    pr = cp * mu / k
    # Smooth-channel reference correlation; corrugated chevron correlations need measured geometry/calibration.
    nu, regime = channel_nusselt(dp.reynolds, pr)
    h_secondary = nu * k / dh
    area = geometry.active_area_m2
    r_total = 1/(refrigerant_h_w_m2_k*area) + geometry.wall_thickness_m/(geometry.wall_conductivity_w_m_k*area) + 1/(h_secondary*area) + fouling_resistance_m2k_w/area
    ua = 1/r_total
    screen = constant_saturation_evaporator(ua_w_k=ua, capacity_rate_w_k=secondary_mass_flow_kg_s*cp,
        secondary_inlet_c=secondary_inlet_c, saturation_temperature_c=refrigerant_saturation_c)
    return PlateEvaporatorResult(h_secondary, dp.reynolds, nu, dp.pressure_drop_pa, regime, ua, screen,
        f"{regime.capitalize()} flow. Preliminary smooth-channel secondary-fluid correlation; refrigerant boiling h is prescribed. Real chevron plates reach turbulence near Re ~ 100-400 and give several times higher h, so this is conservative. Do not use for final selection without a documented chevron/boiling correlation and validation.")
