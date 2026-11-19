# Verification, validation, and limits

This page separates three kinds of evidence:

1. **Verification:** the equations are implemented correctly.
2. **Independent data check:** the thermodynamics agree with published manufacturer data.
3. **Hardware validation:** system predictions agree with measured equipment. This is not done yet.

## 1. Verification (automated, `pytest`)

Run `python -m pytest`. The suite has 59 tests, including:

| Check | Reference |
|---|---|
| Cycle first law `Q_c = Q_e + W` | to numerical precision |
| Valve `h4 = h3`, compressor `η_is` recovered from states | definitions |
| Cycle COP below Carnot COP | second law |
| LMTD with terminal ΔT 50 K / 30 K | 39.1523 K (hand calculation) |
| Counterflow ε-NTU, `C_r = 1`, `NTU = 2` | ε = 2/3 (closed form) |
| UA screening at 950 W/K, C = 1000 W/K, 12 → 4 °C | ε = 1 − e^−0.95, Q = 4906.07 W |
| Water cooling task, 0.239 kg/s, 12 → 7 °C | m·cp·ΔT ≈ 5.02 kW |
| Plate channel laminar Nu | 7.54 (parallel plates), not 4.36 (tube) |
| Gnielinski Nu at Re = 10⁴, Pr = 7 | 79.5 (hand calculation) |
| Coupled solver: both UA residuals | < 10⁻⁵, all five refrigerants |
| Secondary-stream energy balances in the coupled solution | exact |
| Infeasible cases (undersized exchangers, near-critical sink) | clear `InfeasibleOperatingPoint` |
| Compressor map rows with COP above Carnot | rejected (catches Btu/h vs W) |
| Map envelope | convex hull, no extrapolation |
| Streamlit app | every tab clicked through headlessly |

These tests confirm that the code solves the stated equations. They do not show that the equations describe a real machine.

## 2. Independent check against a manufacturer datasheet

**Source:** Copeland ZR34K3E-TFD scroll compressor, R407C, 60 Hz. Sheet 35166-460, printed 2025-09-13. Rated at 20 °F superheat, 15 °F subcooling, dew-point temperatures. 57 published points, from −10 to 55 °F evaporating and 80 to 150 °F condensing. The transcription is in `examples/copeland_zr34k3e_tfd_datasheet.csv`. Each row was checked so that capacity ÷ power reproduces the printed EER within rounding.

The sheet publishes **refrigerant mass flow** and **isentropic efficiency** at every point. Neither is an input to this model, so both are independent tests:

| Quantity | How the model predicts it | Mean abs. error | Worst point |
|---|---|---|---|
| Mass flow | `m = C / (h1 − h4)` at rating superheat/subcooling | **0.14 %** | 0.40 % |
| Isentropic efficiency | `η = M (h2s − h1) / P` | **0.09 points** | 0.31 points |

Sensitivity of this check to modelling choices:

| Modelling choice | Mass-flow error |
|---|---|
| Dew-point convention + datasheet rating conditions (current) | 0.14 % |
| Bubble-point condensing pressure (previous version) | 4.6 % |
| Wrong superheat/subcooling (5 K / 3 K instead of 11.1 K / 8.3 K) | 9.0 % |

So the check is sharp enough to catch real convention errors. It found and confirmed the dew-point fix.

**What this does and does not show:** it verifies the CoolProp property model, the state-point construction, and the dew/bubble and superheat/subcooling handling. It does **not** validate the heat-exchanger UA models or the coupled operating-point prediction. Those need measured system data.

### Compressor-map fit (held-out)

The fitted second-order map is scored only on rows it never saw (75/25 split, seed 2026):

| Output | Held-out MAPE | Worst held-out point |
|---|---|---|
| Capacity | 0.72 % | 2.5 % |
| Electrical power | 0.48 % | 1.1 % |

## 3. Hardware validation (not yet done)

The **Test Validation** tab scores predicted against measured capacity and power for independent test points (MAE, RMSE, MAPE, bias, max error). Doing this properly needs measured system data from a test rig, collected separately from any calibration data. See `REAL_DATA_WORKFLOW.md`.

## Numerical convergence

The moving-boundary exchanger is discretised into 80 segments. Doubling to 160 changes the default case's capacity, COP, and saturation temperatures by less than 0.05 %. The condenser minimum approach changes by 0.8 %, because the pinch sits at the desuperheating/condensing boundary, which the grid only resolves to one segment width. Use more segments if a design decision depends on the exact pinch value.

## Known limitations

- **Steady state only:** no start-up, cycling, charge migration, or control dynamics.
- **No refrigerant-side pressure drop:** saturation pressure is constant through each exchanger.
- **Constant UA:** real coefficients vary with flow, quality, and surface condition.
- **Prescribed refrigerant coefficients:** boiling and condensing h are inputs, not predicted.
- **Conservative plate model:** the secondary side uses smooth-channel correlations, and chevron plates transfer more heat.
- **Dry air only:** no dehumidification (wet-coil) model.
- **Subcritical only:** no transcritical CO₂ gas cooler.
- **Single-speed map:** the example map is at one speed (60 Hz), and variable speed needs multi-speed data.

## Appropriate use

Education, portfolio demonstration, preliminary trade studies, and as a foundation for calibrated engineering software. It is not a substitute for a certified rating, hardware selection by the equipment manufacturer, or safety review.
