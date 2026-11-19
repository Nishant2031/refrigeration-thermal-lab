# Step-by-Step Teaching Guide

## Learning goal

By the end, you should be able to trace every energy flow, explain why the
four refrigerant states move where they do on a pressure-enthalpy diagram, and
show how finite heat-exchanger size determines the cycle's actual operating
temperatures.

## Lesson 1 — Read one thermodynamic state

Open `src/thermal_model/state.py`. CoolProp needs two independent properties.
For example, `(P,T)` defines a superheated or subcooled state, while `(P,H)` is
useful after a valve or compressor calculation. `StatePoint` then collects `P`,
`T`, `h`, `s`, density, quality (when meaningful), and phase.

Exercise: run the example, open `results/cycle_state_points.csv`, and identify
which states are vapor, liquid, and two-phase. Notice that quality is meaningful
only inside the saturation dome.

## Lesson 2 — Build the four-state cycle by hand

Open `src/thermal_model/cycle.py` and follow this order:

1. Convert saturation temperatures to low/high pressures.
2. Add superheat to define state 1.
3. Find the ideal isentropic state 2s at `(p_high,s_1)`.
4. Correct the enthalpy rise using compressor efficiency to obtain state 2.
5. Remove subcooling from the condensing temperature to define state 3.
6. Hold enthalpy constant through the valve to define state 4.
7. Apply the three component energy balances.

Checkpoint: explain why a lower compressor efficiency changes state 2 but does
not directly change states 1, 3, and 4 when saturation temperatures are fixed.

## Lesson 3 — Verify the first law

For the refrigerant control volume, calculate `Q_c - Q_e - W_c`. It should be
zero to floating-point precision. Then confirm `COP_heating = COP_cooling + 1`.
These are identities, not experimental validation; they prove the bookkeeping is
internally consistent.

Exercise: change mass flow from 0.03 to 0.06 kg/s. Duties and power should double,
but COP and pressure ratio should remain unchanged in the standalone cycle.

## Lesson 4 — Understand LMTD

LMTD is the correct mean driving temperature for an exchanger whose local duty
obeys `dQ = U DeltaT dA`. Open `heat_exchangers.py::lmtd`.

Exercise: evaluate `lmtd(50,30)`. The answer is about 39.152 K, between the two
terminal differences. Try a negative terminal difference: the model rejects the
temperature cross instead of returning a misleading number.

## Lesson 5 — Understand effectiveness-NTU

Use ε-NTU when UA and inlet conditions are known but outlet temperatures are not.
`Q_max` is the most the smaller-capacity stream could exchange without violating
the inlet-temperature limit. Effectiveness is the fraction actually achieved.

Exercise: hold both capacity rates fixed and increase UA. Watch effectiveness and
duty rise asymptotically. Compare counterflow and parallel flow at equal NTU.

## Lesson 6 — Phase change is a special capacity-rate limit

During ideal boiling or condensing at fixed pressure, refrigerant temperature is
constant while it absorbs/releases latent heat. Mathematically its capacity rate
tends to infinity and `C_r` tends to zero. This produces the simple relation
`epsilon=1-exp(-NTU)` used by `phase_change_hx`.

Exercise: in an evaporator, increase source inlet temperature while keeping UA,
source capacity rate, and saturation temperature fixed. Duty grows linearly with
the inlet temperature difference.

## Lesson 7 — Add the real refrigerant zones

An actual condenser first desuperheats vapor, then condenses it, then subcools
liquid. An evaporator boils the mixture and usually superheats vapor. Open
`moving_boundary.py`: its pressure/enthalpy nodes let CoolProp locate all those
zones without manually guessing their boundaries.

Open `results/coupled_hx_temperature_profiles.png`. Flat refrigerant portions are
two-phase regions; sloped portions are single-phase. The narrowest vertical gap
is the minimum approach temperature (pinch).

Exercise: compare results using 20, 40, 80, and 160 segments. Required UA should
converge. Eighty is a good teaching/default compromise, not a universal optimum.

## Lesson 8 — Couple the cycle and exchangers

A standalone cycle can accept arbitrary saturation temperatures. Real equipment
cannot: finite UA and external-fluid temperatures decide how far apart the fluids
must be. `coupled.py` guesses `T_e` and `T_c`, solves the entire cycle, computes the
required UA of both exchangers, and iterates until required equals available.

Exercise A: increase condenser UA. Condensing temperature and pressure ratio
should fall, while COP should rise.

Exercise B: decrease source capacity rate. The source experiences a larger
temperature change; eventually a pinch may make the requested duty impossible.

Exercise C: lower compressor efficiency. Compressor outlet temperature and power
rise. The condenser must reject the extra work, so the coupled operating point
also shifts.

## Lesson 9 — Read the plots like an engineer

- `cycle_ph_diagram.png`: valve line is constant enthalpy; evaporation and
  condensation occur mostly inside the dome.
- `cop_vs_condensing_temperature.png`: increasing lift costs compressor work and
  reduces COP; the Carnot curve is an unattainable upper bound here.
- `coupled_hx_temperature_profiles.png`: check approach temperatures and zones;
  a visually tiny gap signals large required area or a pinch.

## Lesson 10 — Extend responsibly

Recommended next additions, in order:

1. heat-exchanger and line pressure drops;
2. compressor map efficiency and volumetric displacement;
3. temperature-dependent secondary-fluid heat capacity;
4. correlations for refrigerant-side and secondary-side heat-transfer coefficients;
5. charge/inventory and transient wall-fluid energy storage;
6. calibration against manufacturer or laboratory data with uncertainty bands.

For each addition, write one conservation test, one analytical limiting-case test,
and one expected-trend test before trusting new plots.

