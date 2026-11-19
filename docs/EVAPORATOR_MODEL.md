# Evaporator module

## What is implemented

The **3 · Evaporator** tab intentionally keeps three questions separate:

1. **Cooling requirement** - calculates the duty a liquid stream needs to reach a requested outlet temperature. With a declared constant specific heat,
   `Q_required = m_dot * cp * (T_in - T_target)`.
2. **Constant-saturation UA screening** - calculates the duty a finite-UA exchanger can deliver when refrigerant saturation temperature is prescribed:
   `epsilon = 1 - exp(-UA / C_secondary)`,
   `Q = epsilon * C_secondary * (T_secondary,in - T_sat)`.
3. **Refrigerant-state check** - uses CoolProp to calculate saturation pressure, inlet enthalpy from quality, outlet enthalpy from pressure/temperature, and
   `Q_ref = m_dot_ref * (h_out - h_in)` for pure refrigerants. An outlet below saturation is rejected (impossible without pressure drop); an outlet at saturation is treated as saturated vapor; zeotropic blends are detected from their glide and rejected.

The first is a load requirement; the second is a bounded heat-transfer screening result; the third is a property and energy-balance check. None independently rates an arbitrary physical evaporator.

## Verification fixtures

Automated tests cover water cooling loads (2, 5, 7 kW), the 950 W/K constant-saturation screening case (4906.07 W, checked against the closed-form ε-NTU result), and the R134a property case at 5 °C saturation. Run all checks with `python -m pytest`.

## Not yet implemented

The application does **not** yet claim a high-fidelity physical rating for a named plate or finned-tube evaporator. That needs a selected geometry, a documented two-phase boiling heat-transfer correlation, compatible friction/acceleration pressure-drop correlation, refrigerant distribution model, local segment iteration, correlation validity checks, and calibration/independent validation data. These are the next engineering steps for a geometry-specific rating model.
