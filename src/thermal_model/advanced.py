"""Engineering-estimation extensions with explicit applicability limits.

These routines support preliminary/design-model estimates. They are deliberately
separate from the validated thermodynamic core so assumptions remain visible.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, log10, pi


@dataclass(frozen=True)
class TubeBundleGeometry:
    """Simple smooth-tube exchanger geometry; SI units."""

    tube_inner_diameter_m: float
    tube_outer_diameter_m: float
    tube_length_m: float
    tube_count: int
    wall_conductivity_w_m_k: float = 385.0  # copper
    fouling_resistance_inner_m2k_w: float = 0.0
    fouling_resistance_outer_m2k_w: float = 0.0

    @property
    def inner_area_m2(self) -> float:
        return pi * self.tube_inner_diameter_m * self.tube_length_m * self.tube_count

    @property
    def outer_area_m2(self) -> float:
        return pi * self.tube_outer_diameter_m * self.tube_length_m * self.tube_count


@dataclass(frozen=True)
class UAEstimate:
    ua_w_k: float
    overall_u_outer_w_m2_k: float
    inner_area_m2: float
    outer_area_m2: float
    resistance_inner_k_w: float
    resistance_wall_k_w: float
    resistance_outer_k_w: float
    resistance_fouling_k_w: float


def estimate_tube_bundle_ua(
    geometry: TubeBundleGeometry, h_inner_w_m2_k: float, h_outer_w_m2_k: float,
) -> UAEstimate:
    """Overall-UA from cylindrical-wall thermal resistances.

    Heat-transfer coefficients must originate from an applicable, documented
    correlation or measured data. Do not use this routine to certify a coil.
    """
    if min(geometry.tube_inner_diameter_m, geometry.tube_outer_diameter_m,
           geometry.tube_length_m, geometry.tube_count, geometry.wall_conductivity_w_m_k,
           h_inner_w_m2_k, h_outer_w_m2_k) <= 0:
        raise ValueError("Geometry, conductivity, and film coefficients must be positive.")
    if geometry.tube_outer_diameter_m <= geometry.tube_inner_diameter_m:
        raise ValueError("Outer diameter must exceed inner diameter.")
    ai, ao = geometry.inner_area_m2, geometry.outer_area_m2
    r_i = 1.0 / (h_inner_w_m2_k * ai)
    # Cylindrical wall conduction: ln(Do/Di) / (2*pi*k*L_total)
    r_w = log(geometry.tube_outer_diameter_m / geometry.tube_inner_diameter_m) / (2 * pi * geometry.wall_conductivity_w_m_k * geometry.tube_length_m * geometry.tube_count)
    r_o = 1.0 / (h_outer_w_m2_k * ao)
    r_f = geometry.fouling_resistance_inner_m2k_w / ai + geometry.fouling_resistance_outer_m2k_w / ao
    r_total = r_i + r_w + r_o + r_f
    ua = 1.0 / r_total
    return UAEstimate(ua, ua / ao, ai, ao, r_i, r_w, r_o, r_f)


@dataclass(frozen=True)
class PressureDropResult:
    pressure_drop_pa: float
    reynolds: float
    friction_factor: float
    velocity_m_s: float
    regime: str


def single_phase_pressure_drop(
    *, mass_flow_kg_s: float, density_kg_m3: float, viscosity_pa_s: float,
    hydraulic_diameter_m: float, flow_area_m2: float, length_m: float,
    roughness_m: float = 1.5e-6, minor_loss_coefficient: float = 0.0,
    laminar_friction_constant: float = 64.0,
) -> PressureDropResult:
    """Darcy-Weisbach pressure drop for a single-phase pipe/channel.

    ``laminar_friction_constant`` is the fully developed laminar Darcy
    constant C in f = C/Re: 64 for a round tube, 96 for flow between wide
    parallel plates (Shah & London).  Two-phase pressure drop requires a
    regime-specific correlation and is not inferred from this function.
    """
    if min(mass_flow_kg_s, density_kg_m3, viscosity_pa_s, hydraulic_diameter_m,
           flow_area_m2, length_m, laminar_friction_constant) <= 0 or roughness_m < 0 or minor_loss_coefficient < 0:
        raise ValueError("Inputs must be physically positive (roughness/loss can be zero).")
    velocity = mass_flow_kg_s / (density_kg_m3 * flow_area_m2)
    re = density_kg_m3 * velocity * hydraulic_diameter_m / viscosity_pa_s
    if re < 2300:
        f = laminar_friction_constant / re
        regime = "laminar"
    else:
        # Swamee-Jain explicit Darcy-friction estimate for turbulent internal flow.
        f = 0.25 / (log10(roughness_m / (3.7 * hydraulic_diameter_m) + 5.74 / re**0.9) ** 2)
        regime = "turbulent"
    dynamic = 0.5 * density_kg_m3 * velocity**2
    dp = (f * length_m / hydraulic_diameter_m + minor_loss_coefficient) * dynamic
    return PressureDropResult(dp, re, f, velocity, regime)


@dataclass(frozen=True)
class CompressorEstimate:
    refrigerant_mass_flow_kg_s: float
    volumetric_efficiency: float
    isentropic_efficiency: float
    shaft_power_w: float
    electrical_power_w: float


def compressor_displacement_estimate(
    *, suction_density_kg_m3: float, displacement_m3_rev: float, speed_rpm: float,
    pressure_ratio: float, isentropic_specific_work_j_kg: float,
    motor_efficiency: float = 0.92, volumetric_efficiency_at_ratio_1: float = 0.90,
    volumetric_efficiency_slope: float = 0.035, isentropic_efficiency: float = 0.72,
) -> CompressorEstimate:
    """Transparent displacement-compressor estimate, not a manufacturer map."""
    if min(suction_density_kg_m3, displacement_m3_rev, speed_rpm, pressure_ratio,
           isentropic_specific_work_j_kg, motor_efficiency, isentropic_efficiency) <= 0:
        raise ValueError("Compressor inputs must be positive.")
    eta_v = max(0.35, min(1.0, volumetric_efficiency_at_ratio_1 - volumetric_efficiency_slope * (pressure_ratio - 1)))
    mdot = suction_density_kg_m3 * displacement_m3_rev * speed_rpm / 60.0 * eta_v
    shaft = mdot * isentropic_specific_work_j_kg / isentropic_efficiency
    return CompressorEstimate(mdot, eta_v, isentropic_efficiency, shaft, shaft / motor_efficiency)


@dataclass(frozen=True)
class AuxiliaryPowerResult:
    pump_power_w: float
    fan_power_w: float
    total_auxiliary_power_w: float


def auxiliary_power(*, liquid_pressure_drop_pa: float = 0.0, liquid_volume_flow_m3_s: float = 0.0,
                    pump_efficiency: float = 0.65, air_pressure_drop_pa: float = 0.0,
                    air_volume_flow_m3_s: float = 0.0, fan_efficiency: float = 0.55) -> AuxiliaryPowerResult:
    if pump_efficiency <= 0 or fan_efficiency <= 0:
        raise ValueError("Pump and fan efficiency must be positive.")
    pump = liquid_pressure_drop_pa * liquid_volume_flow_m3_s / pump_efficiency
    fan = air_pressure_drop_pa * air_volume_flow_m3_s / fan_efficiency
    if pump < 0 or fan < 0:
        raise ValueError("Pressure losses and volume flows cannot be negative.")
    return AuxiliaryPowerResult(pump, fan, pump + fan)


@dataclass(frozen=True)
class ChargeEstimate:
    refrigerant_charge_kg: float
    liquid_fraction_mass: float
    vapor_fraction_mass: float


def homogeneous_two_phase_charge(*, volume_m3: float, quality: float, liquid_density_kg_m3: float,
                                 vapor_density_kg_m3: float) -> ChargeEstimate:
    """Homogeneous-equilibrium inventory estimate for a two-phase control volume."""
    if volume_m3 <= 0 or liquid_density_kg_m3 <= 0 or vapor_density_kg_m3 <= 0 or not 0 <= quality <= 1:
        raise ValueError("Use positive volume/densities and quality in [0, 1].")
    mixture_density = 1.0 / ((1 - quality) / liquid_density_kg_m3 + quality / vapor_density_kg_m3)
    mass = mixture_density * volume_m3
    return ChargeEstimate(mass, 1 - quality, quality)

