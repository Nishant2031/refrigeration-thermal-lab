# Real Hardware Data Workflow

## 1. Compressor-map CSV

Required columns, with SI-based units shown in names:

```csv
evaporating_temperature_c,condensing_temperature_c,speed_rpm,capacity_w,power_w
```

Each row must describe one distinct, stable operating point for the **same**
compressor model, refrigerant, oil, electrical configuration, and test method.
Do not combine data from different compressor variants or refrigerants.

Example:

```csv
evaporating_temperature_c,condensing_temperature_c,speed_rpm,capacity_w,power_w
-10,35,3600,3900,1100
-5,40,3600,4500,1280
0,45,3600,4900,1460
```

The app validates the file, randomly but reproducibly reserves a hold-out set,
fits a transparent second-order response surface to the remaining map points,
and reports error only on the hold-out points. It requires at least 14 distinct
map rows. Use substantially more for a credible map.

The fitted response surface is a convenience representation of the imported
map, not a new manufacturer rating. The app refuses to use it outside the map's
envelope, defined as the convex hull of the published (evaporating, condensing)
points within the speed range. Real envelopes are not rectangles.

**Units matter.** Capacity and power must both be in W. US datasheets publish
capacity in Btu/h (1 Btu/h = 0.29307107 W). Every row is checked against the
Carnot limit `COP < T_e / (T_c - T_e)`, which rejects a Btu/h-labelled-as-W
map immediately.

**Rating conditions.** Record the map's rating superheat and subcooling (for
example 20 °F / 15 °F = 11.1 K / 8.3 K). The coupled solver converts map
capacity to mass flow at those conditions, then corrects to your operating
superheat and subcooling.

The bundled example, `examples/copeland_zr34k3e_tfd_map.csv`, is the complete
57-point Copeland ZR34K3E-TFD R407C sheet. `copeland_zr34k3e_tfd_datasheet.csv`
keeps every published column in original units, including mass flow and
isentropic efficiency, which the datasheet check in tab 6 uses.

## 2. Evaporator and condenser data

Enter or retain these facts for each coil:

- manufacturer, model, serial number, drawing revision, and photograph;
- refrigerant, circuit count, tube inner/outer diameter, tube length, tube count;
- fin geometry or water-side geometry, material, and fouling condition;
- airflow/water-flow data and pressure drop;
- installed sensors, locations, calibration certificate IDs, and acquisition rate;
- selected heat-transfer / pressure-drop correlation and its published validity range.

The **4 · Coil Geometry** tab calculates a transparent cylindrical-wall UA resistance
network. Its refrigerant-side coefficient is an input until a correlation has
been selected and implemented for the exact coil and flow regime.

## 3. Test-point CSV

After locking a model configuration, export its predictions and compare them
with measurements from **independent** physical tests:

```csv
test_id,predicted_capacity_w,measured_capacity_w,predicted_power_w,measured_power_w
TP-01,4500,4380,1280,1310
```

The interface calculates MAE, RMSE, MAPE, bias, and maximum absolute percent
error separately for capacity and electrical power.

Do not validate using the same points used for calibration. Keep at least one
test envelope separate: different combinations of evaporating condition,
condensing condition, speed, and flow rate.

## 4. Recommended collection order

1. Obtain the compressor manufacturer map and its test-condition documentation.
2. Document exact coil geometry and refrigerant circuit arrangement.
3. Measure steady-state baseline points with calibrated instruments.
4. Fit uncertain parameters using only the designated calibration subset.
5. Freeze the parameter set and software commit.
6. Run independent validation points and report errors and uncertainty.
7. Only then define whether the model is acceptable inside its demonstrated
   operating envelope.

