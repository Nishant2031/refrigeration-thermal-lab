import pytest

from thermal_model.advanced import (
    TubeBundleGeometry, auxiliary_power, compressor_displacement_estimate,
    estimate_tube_bundle_ua, homogeneous_two_phase_charge, single_phase_pressure_drop,
)


def test_tube_bundle_ua_has_positive_resistances_and_area():
    result = estimate_tube_bundle_ua(TubeBundleGeometry(0.008, 0.010, 4.0, 24), 1600, 120)
    assert result.ua_w_k > 0
    assert result.inner_area_m2 > 0
    assert result.outer_area_m2 > result.inner_area_m2
    assert result.resistance_inner_k_w + result.resistance_wall_k_w + result.resistance_outer_k_w + result.resistance_fouling_k_w == pytest.approx(1 / result.ua_w_k)


def test_single_phase_pressure_drop_increases_with_flow():
    low = single_phase_pressure_drop(mass_flow_kg_s=0.02, density_kg_m3=998, viscosity_pa_s=0.001, hydraulic_diameter_m=0.01, flow_area_m2=7.85e-5, length_m=8)
    high = single_phase_pressure_drop(mass_flow_kg_s=0.04, density_kg_m3=998, viscosity_pa_s=0.001, hydraulic_diameter_m=0.01, flow_area_m2=7.85e-5, length_m=8)
    assert high.pressure_drop_pa > low.pressure_drop_pa > 0
    assert high.reynolds > low.reynolds


def test_compressor_auxiliary_and_charge_estimates_are_physical():
    c = compressor_displacement_estimate(suction_density_kg_m3=18, displacement_m3_rev=12e-6, speed_rpm=3600, pressure_ratio=3, isentropic_specific_work_j_kg=25_000)
    a = auxiliary_power(liquid_pressure_drop_pa=30_000, liquid_volume_flow_m3_s=0.0007, air_pressure_drop_pa=160, air_volume_flow_m3_s=0.7)
    q = homogeneous_two_phase_charge(volume_m3=0.0015, quality=0.25, liquid_density_kg_m3=1200, vapor_density_kg_m3=25)
    assert c.electrical_power_w > c.shaft_power_w > 0
    assert 0 < c.volumetric_efficiency <= 1
    assert a.total_auxiliary_power_w == pytest.approx(a.pump_power_w + a.fan_power_w)
    assert 0 < q.refrigerant_charge_kg < 1.8

