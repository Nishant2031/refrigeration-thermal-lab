from thermal_model.plate_evaporator import PlateEvaporatorGeometry, plate_evaporator_screening


def test_plate_screening_returns_physical_preliminary_result():
    result = plate_evaporator_screening(geometry=PlateEvaporatorGeometry(.5, .12, 30, .0025),
        secondary_fluid="Water", secondary_mass_flow_kg_s=.24, secondary_inlet_c=12,
        refrigerant_saturation_c=4, refrigerant_h_w_m2_k=1800)
    assert result.ua_w_k > 0
    assert result.screening.duty_w > 0
    assert 4 < result.screening.secondary_outlet_c < 12
    assert result.secondary_pressure_drop_pa > 0
