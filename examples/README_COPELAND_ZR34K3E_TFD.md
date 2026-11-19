# Copeland ZR34K3E-TFD: how to import this datasheet

## What was read from the supplied sheet

| Item | Value |
|---|---|
| Manufacturer | Copeland |
| Exact model | ZR34K3E-TFD |
| Mechanism | Scroll |
| Refrigerant | HFC-407C |
| Map basis | Saturated dew-point evaporating and condensing temperatures |
| Frequency | 60 Hz |
| Electrical configuration | 460 V, three phase |
| Rating superheat | 20 F (11.1 K) |
| Rating subcooling | 15 F (8.3 K) |
| Condenser ambient | 95 F (35 C) |
| Sheet accuracy statement | Nominal values, plus or minus 5 percent, based on 72-hour run-in |

The sheet columns mean: `C` capacity in Btu/h, `P` electrical power in W, `A` current in A, `M` refrigerant mass flow in lb/h, `E` EER in Btu/Wh, and `%` isentropic efficiency.

## Use it in the app

1. In **Compressor Map**, click **Load bundled Copeland R407C example** (or upload this CSV).
2. In **Coupled Solver**, select `R407C` and set map rating superheat 11.1 K / subcooling 8.3 K. Do not use this map for R134a, R32, R290, or another refrigerant.
3. Keep the metadata JSON with the CSV. The map values only make sense with the stated 20 F superheat, 15 F subcooling, 95 F ambient, and 60 Hz conditions.
4. Fit the map and inspect held-out errors. It is an interpolation check, not yet a complete-system validation.

## Important limits

- `copeland_zr34k3e_tfd_datasheet.csv` contains all 57 published points (condensing 80-150 F, evaporating -10 to 55 F) with every column: capacity, power, current, mass flow, EER, and isentropic efficiency, in original units. Each row was checked so that capacity/power reproduces the printed EER within rounding.
- `copeland_zr34k3e_tfd_map.csv` is the same 57 points in the app's SI map schema. The envelope is not rectangular (high condensing temperatures are only published at warm evaporating temperatures); the app uses the convex hull of these points as the envelope.
- Temperatures are converted to Celsius. **Capacity is converted from the sheet's Btu/h to W** (`capacity_w = capacity_btu_h × 0.29307107`); the original Btu/h value is kept in the `capacity_btu_h` column for traceability. Power is in W as published.
- An earlier revision of this CSV stored the Btu/h numbers in the `capacity_w` column. That implied COPs above the Carnot limit; the app now rejects any map row whose COP exceeds Carnot.
- `3600 rpm` is a nominal representation of the sheet's 60 Hz operation because the current app expects speed in rpm. It is not a documented shaft-speed measurement.
- The map must not be extrapolated outside these imported conditions. Obtain the full manufacturer export and operating-limit chart before selecting a real system.
