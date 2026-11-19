# Equations, Assumptions, and Numerical Method

## 1. State convention

The loop is `1 -> 2 -> 3 -> 4 -> 1`:

1. compressor inlet: low-pressure superheated vapor;
2. compressor outlet: high-pressure superheated vapor;
3. condenser outlet: high-pressure subcooled liquid;
4. valve outlet: low-pressure two-phase mixture.

Evaporating and condensing inputs are refrigerant saturation temperatures. The
saturation pressures are obtained from CoolProp:

```text
p_e = p_sat(T_e),   p_c = p_sat(T_c)
T_1 = T_e + superheat
T_3 = T_c - subcooling
```

## 2. Compressor and valve

The isentropic reference outlet is found at `(p_c, s_1)`:

```text
h_2s = h(p_c, s_1)
eta_is = (h_2s - h_1) / (h_2 - h_1)
h_2 = h_1 + (h_2s - h_1) / eta_is
```

The expansion valve is adiabatic, has no shaft work, and is isenthalpic:

```text
h_4 = h_3
```

## 3. Cycle energy balances

For steady refrigerant mass flow `m_dot`:

```text
Q_e = m_dot (h_1 - h_4)
W_c = m_dot (h_2 - h_1)
Q_c = m_dot (h_2 - h_3)
Q_c = Q_e + W_c
COP_cooling = Q_e / W_c
COP_heating = Q_c / W_c = COP_cooling + 1
pressure ratio = p_c / p_e
```

The ideal reversed-Carnot upper bound between saturation temperatures is:

```text
COP_Carnot = T_e / (T_c - T_e)       [absolute temperatures]
```

## 4. LMTD

For positive terminal temperature differences:

```text
DeltaT_lm = (DeltaT_1 - DeltaT_2) / ln(DeltaT_1 / DeltaT_2)
Q = U A F DeltaT_lm
```

`F=1` for the modeled one-pass counterflow and parallel-flow arrangements. The
code uses the exact limiting value when both terminal differences are equal and
rejects a temperature cross. A shell-and-tube correction factor is not silently
assumed; add an explicit `F` model for multipass equipment.

## 5. Sensible/sensible effectiveness-NTU

```text
C_h = m_dot_h cp_h                  C_c = m_dot_c cp_c
C_min = min(C_h, C_c)               C_r = C_min / C_max
NTU = UA / C_min
Q_max = C_min (T_h,in - T_c,in)     Q = epsilon Q_max
```

Counterflow:

```text
epsilon = [1 - exp(-NTU(1-C_r))] / [1 - C_r exp(-NTU(1-C_r))]
```

At `C_r=1`, its analytical limit is `epsilon=NTU/(1+NTU)`.

Parallel flow:

```text
epsilon = [1 - exp(-NTU(1+C_r))] / (1+C_r)
```

The code also supplies the standard both-unmixed crossflow correlation.

## 6. Isothermal phase-change exchanger

For a refrigerant fixed at `T_sat`, its effective capacity rate is infinite, so
`C_r -> 0`. The sensible stream is `C_min`:

```text
NTU = UA / C_sensible
epsilon = 1 - exp(-NTU)
Q = epsilon C_sensible |T_sensible,in - T_sat|
```

This is appropriate for a phase-change-only zone. It does not include
superheating, desuperheating, or subcooling; the moving-boundary model does.

## 7. Moving-boundary exchanger

The coupled model divides the refrigerant enthalpy change into `N` equal-duty
segments. At each node CoolProp evaluates:

```text
T_r,i = T(p, h_i)
```

This automatically yields a nearly constant temperature in the two-phase zone
and changing temperature in single-phase zones. The counterflow secondary-fluid
temperature follows its energy balance. For each segment:

```text
DeltaQ_i = m_dot_r |h_(i+1)-h_i|
UA_i = DeltaQ_i / DeltaT_lm,i
UA_required = sum(UA_i)
```

If any local temperature difference is nonpositive, the requested duty is
thermally impossible for those inlet conditions and the model reports a pinch.

## 8. Coupled solution

The two unknowns are `T_e` and `T_c`. SciPy's bounded least-squares solver drives:

```text
r_e = (UA_e,required - UA_e,available) / UA_e,available
r_c = (UA_c,required - UA_c,available) / UA_c,available
```

to zero. Every solver evaluation recalculates the complete cycle and both heat
exchangers. Thus changing compressor efficiency, mass flow, inlet temperatures,
capacity rates, or either UA changes the operating saturation temperatures.

**Robust closure (`solver.py`).** The least-squares solve is tried from three
starting points. If none converges, a bracketed nested root search is used. It
relies on two physical monotonicities: the condenser's required UA falls as `T_c`
rises, and the evaporator's required UA rises as `T_e` rises. A pinch is treated
as "infinite UA required", which keeps both functions monotone. If no bracket
exists inside the allowed window, the case is physically infeasible and the solver
raises `InfeasibleOperatingPoint`, naming the limiting exchanger.

The allowed window is `T_e` in [max(T_source - 40, T_min + 5), T_source - 0.2] and
`T_c` in [T_sink + 0.2, min(T_sink + 50, T_crit - 3)]. For compressor-map closure
it is further intersected with the map envelope, so the map is never extrapolated.

**Map closure correction.** Manufacturer capacity is converted to mass flow
with the refrigerating effect at the map's *rating* superheat and subcooling, then
the actual capacity uses the *operating* superheat and subcooling (AHRI 540 practice):

```text
m_dot = Q_map / (h1 - h4)_rating
Q_actual = m_dot (h1 - h4)_operating
```

## 8a. Saturation-temperature convention

Evaporating and condensing temperatures are dew-point temperatures, matching
compressor-rating standards (AHRI 540, EN 12900). For pure fluids, dew and bubble
points coincide. For zeotropic blends such as R407C, the high-side pressure is the
dew pressure at `T_c`, and subcooling is measured from the bubble temperature at
that pressure.

## 8b. Plate-channel secondary side

Hydraulic diameter `D_h = 2b` for a wide channel of gap `b`.

```text
Laminar  (Re < 2300):  Nu = 7.54          (parallel plates, both walls at uniform T)
                       f  = 96 / Re       (Darcy, parallel plates)
Turbulent (Re >= 1e4): Nu = Gnielinski = (f/8)(Re-1000)Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3)-1)),
                       f  = (0.790 ln Re - 1.64)^-2
Transition:            linear interpolation in Re between Nu(2300) and Nu(1e4)
```

The round-tube values (Nu = 4.36, f = 64/Re) do not apply to a plate channel.
Real chevron plates become turbulent at much lower Re and transfer more heat, so
the smooth-channel result is a conservative screening value.

## 8c. Data checks

Every compressor-map row must satisfy `COP = Q/W < T_e/(T_c - T_e)` (Carnot).
A violation almost always means capacity in Btu/h labelled as W.

The map envelope is the convex hull of the published (T_e, T_c) points. A fitted
map is never evaluated outside it.

## 8d. Datasheet check

With rating superheat `SH` and subcooling `SC` from the sheet:

```text
m_model   = C_sheet / (h1(T_e,dew + SH) - h3(T_bubble(p_c) - SC))
eta_model = M_sheet (h2s - h1) / P_sheet
```

These are compared with the sheet's published mass flow and isentropic efficiency (see VALIDATION.md).

## 9. Physical assumptions

- steady state; one-dimensional mean flow;
- thermodynamic equilibrium properties (blends use dew/bubble conventions; see 8a);
- single-stage compressor with constant isentropic efficiency;
- negligible kinetic/potential energy and ambient heat leaks;
- no suction/liquid-line heat exchanger;
- no pressure drop in piping or heat exchangers;
- fixed superheat, subcooling, secondary heat-capacity rates, and total UA;
- counterflow coupled heat exchangers;
- no frost, oil, maldistribution, fan/pump power, or transient inventory effects.

