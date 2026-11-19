"""Refrigeration Thermal Lab - CoolProp engineering workspace.

Run:  streamlit run app_streamlit.py
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

# Make ``thermal_model`` importable without ``pip install -e .`` (e.g. on Streamlit Community Cloud).
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from thermal_model.advanced import TubeBundleGeometry, estimate_tube_bundle_ua
from thermal_model.coupled import CoupledInputs, solve_coupled_system
from thermal_model.data_workflow import (COMPRESSOR_MAP_COLUMNS, TEST_POINT_COLUMNS, ErrorMetrics,
                                         fit_compressor_map, score_test_points, validate_compressor_map)
from thermal_model.datasheet_validation import validate_against_datasheet
from thermal_model.evaporator import (constant_saturation_evaporator, liquid_cooling_task,
                                      refrigerant_evaporator_states)
from thermal_model.map_coupled import MapCoupledInputs, solve_map_coupled_system
from thermal_model.plate_evaporator import PlateEvaporatorGeometry, plate_evaporator_screening
from thermal_model.plots import cycle_ph_figure, hx_profiles_figure
from thermal_model.sweep import SWEEPABLE, diminishing_returns, sweep_coupled

# Streamlit >= 1.50 uses width="stretch"; older versions use use_container_width=True.
_ST_VERSION = tuple(int(x) for x in st.__version__.split(".")[:2] if x.isdigit())
STRETCH = {"width": "stretch"} if _ST_VERSION >= (1, 50) else {"use_container_width": True}

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_MAP = APP_DIR / "examples" / "copeland_zr34k3e_tfd_map.csv"
EXAMPLE_DATASHEET = APP_DIR / "examples" / "copeland_zr34k3e_tfd_datasheet.csv"
EXAMPLE_MAP_INFO = {"name": "Copeland ZR34K3E-TFD scroll datasheet (57 points)", "refrigerant": "R407C",
                    "rating_superheat_k": 11.11, "rating_subcooling_k": 8.33}
REFRIGERANTS = ["R134a", "R407C", "R290", "R32", "R1234yf"]

st.set_page_config(page_title="Refrigeration Thermal Lab", page_icon="❄️", layout="wide")
st.markdown("""
<style>
  .stApp {background:#f5f9fc; color:#102a3a}
  [data-testid="stHeader"] {background:#f5f9fc}
  .block-container {padding-top:2rem; max-width:1380px}
  p, label, li, [data-testid="stMarkdownContainer"] {font-size:1.03rem; color:#17374a}
  [data-testid="stWidgetLabel"] p, [data-testid="stWidgetLabel"] label {font-weight:700!important; color:#102a3a!important; font-size:1rem!important}
  input, textarea, [data-baseweb="select"] > div {background:#ffffff!important; color:#102a3a!important; border-color:#416981!important; font-size:1.03rem!important}
  .hero-label {font: 700 13px monospace; color:#16646d; letter-spacing:.12em}
  .hero-title {font-size:3.35rem; font-weight:800; letter-spacing:-.06em; line-height:1; margin:.25rem 0 .75rem; color:#102a3a}
  .hero-title span {color:#087f8c}
  .hero-copy {color:#294d61; max-width:780px; font-size:1.13rem; line-height:1.6; font-weight:500}
  .cycle-strip {display:grid; grid-template-columns:1fr auto 1fr auto 1fr auto 1fr; align-items:center; gap:8px; margin:1.8rem 0 2rem}
  .cycle-node {border:2px solid #9bb7c7; border-radius:12px; padding:15px 16px; background:#ffffff; min-height:88px}
  .cycle-node b {display:block; color:#102a3a; font-size:1.05rem}.cycle-node small {display:block; color:#35576b; font-size:.95rem; line-height:1.35; margin-top:4px}
  .cycle-arrow {color:#087f8c; font-size:1.6rem; font-weight:800}
  .model-note {border-left:5px solid #087f8c; padding:.85rem 1rem; background:#e7f3f8; color:#16384b; margin:.7rem 0 1.4rem; font-size:1.05rem; border-radius:4px}
  [data-testid="stMetric"] {background:#ffffff; border:2px solid #a3becd; padding:14px; border-radius:10px}
  [data-testid="stMetricLabel"] {color:#17374a!important; font-weight:700!important}
  div[data-testid="stTabs"] button {font-weight:750; color:#23485b; font-size:1rem}
  div[data-testid="stTabs"] button[aria-selected="true"] {color:#075e69}
  [data-testid="stButton"] button {background:#087f8c; color:#ffffff; border:2px solid #05636d; font-size:1.03rem; font-weight:750; min-height:46px}
  [data-testid="stButton"] button:hover, [data-testid="stButton"] button:focus, [data-testid="stButton"] button:active,
  [data-testid="stButton"] button:focus:not(:active) {background:#05636d; color:#ffffff!important; border-color:#05636d}
  [data-testid="stButton"] button *, [data-testid="stButton"] button [data-testid="stMarkdownContainer"] p {color:#ffffff!important}
  @media (max-width:760px) {.hero-title{font-size:2.45rem}.cycle-strip{grid-template-columns:1fr 1fr}.cycle-arrow{display:none}.cycle-node{min-height:96px}}
</style>
<div class="hero-label">VAPOR-COMPRESSION CYCLE · HEAT-EXCHANGER MODEL · COOLPROP</div>
<div class="hero-title">Refrigeration <span>Thermal Lab</span></div>
<div class="hero-copy">Predict where a refrigeration system actually operates from its heat-exchanger sizes and source/sink conditions, size and check the evaporator, and calibrate against real compressor data. Every calculation states its model boundary, so screening estimates are never confused with verified hardware performance.</div>
<div class="cycle-strip">
  <div class="cycle-node"><b>01 · Compressor</b><small>Raises pressure and drives flow</small></div><div class="cycle-arrow">→</div>
  <div class="cycle-node"><b>02 · Condenser</b><small>Rejects heat to the sink</small></div><div class="cycle-arrow">→</div>
  <div class="cycle-node"><b>03 · Expansion valve</b><small>Creates the low-pressure side</small></div><div class="cycle-arrow">→</div>
  <div class="cycle-node"><b>04 · Evaporator</b><small>Absorbs heat from the secondary fluid</small></div>
</div>
<div class="model-note"><b>How to use this app:</b> <b>1 · Coupled Solver</b> finds the steady operating point (evaporating/condensing temperatures, capacity, COP) from exchanger UA values and source/sink conditions.
<b>2 · Design Sweep</b> shows how the system responds when one design variable changes. <b>3 · Evaporator</b> and <b>4 · Coil Geometry</b> size or check the evaporator. <b>5 · Compressor Map</b> and <b>6 · Test Validation</b> bring in real data.
Results are engineering estimates, not a certified rating.</div>
""", unsafe_allow_html=True)


# ----------------------------------------------------------------- helpers
def download_csv(columns: tuple[str, ...], rows: list[dict], label: str):
    st.download_button(label, pd.DataFrame(rows, columns=columns).to_csv(index=False).encode(),
                       label.lower().replace(" ", "_") + ".csv", "text/csv")


def metric_table(metrics: ErrorMetrics) -> pd.DataFrame:
    return pd.DataFrame({"Metric": ["MAE", "RMSE", "MAPE", "Bias", "Max absolute % error"],
                         "Value": [metrics.mae, metrics.rmse, metrics.mape_percent, metrics.bias_percent,
                                   metrics.maximum_absolute_percent_error],
                         "Unit": ["W", "W", "%", "%", "%"]})


def show_figure(fig):
    st.pyplot(fig, **STRETCH)
    plt.close(fig)


def stale_note(key: str, current: dict):
    """Warn when inputs changed after a stored result was calculated."""
    if st.session_state.get(key) not in (None, current):
        st.caption("⚠️ Inputs changed since this result was calculated - run again to update.")


def show_operating_point(result, *, map_driven: bool):
    c = result.cycle
    if not result.solver_success:
        st.warning(f"Solver did not fully converge (UA residuals {result.residual_ua_fraction[0]:.2e}, "
                   f"{result.residual_ua_fraction[1]:.2e}). Treat this point with caution.")
    if map_driven:
        st.success("Map-driven closure: the compressor map sets refrigerant flow and electrical power; "
                   "the heat exchangers close the saturation temperatures.")
        capacity, power = result.map_capacity_w, result.map_electrical_power_w
        cycle_cop, system_cop = result.map_cop, result.system_cop
        aux = result.inputs.auxiliary_power_w
    else:
        capacity, power = c.evaporator_heat_w, c.compressor_power_w
        cycle_cop, system_cop = c.cop_cooling, result.system_cop_cooling
        aux = result.inputs.auxiliary_power_w
    a, b, d, e = st.columns(4)
    a.metric("Cooling capacity", f"{capacity/1000:.2f} kW")
    b.metric("Cycle COP", f"{cycle_cop:.2f}", help="Cooling capacity / compressor power.")
    d.metric("Evaporating temp.", f"{c.inputs.evaporating_temperature_c:.2f} °C", help="Dew-point saturation temperature.")
    e.metric("Condensing temp.", f"{c.inputs.condensing_temperature_c:.2f} °C", help="Dew-point saturation temperature.")
    f, g, h = st.columns(3)
    f.metric("System COP", f"{system_cop:.2f}", help="Cooling capacity / (compressor + declared fan, pump and controls power).")
    g.metric("Carnot COP (same temps)", f"{c.carnot_cop_cooling:.2f}", help="Thermodynamic upper limit between the two saturation temperatures.")
    h.metric("Second-law efficiency", f"{cycle_cop / c.carnot_cop_cooling:.0%}", help="Cycle COP / Carnot COP.")
    rows = {
        "Refrigerant mass flow": f"{c.inputs.refrigerant_mass_flow_kg_s:.4f} kg/s",
        "Compressor power": f"{power/1000:.3f} kW",
        "Declared auxiliary power": f"{aux/1000:.3f} kW",
        "Condenser heat rejected": f"{c.condenser_heat_w/1000:.3f} kW",
        "Pressure ratio": f"{c.pressure_ratio:.2f}",
        "Source outlet": f"{result.evaporator.secondary_outlet_c:.2f} °C",
        "Sink outlet": f"{result.condenser.secondary_outlet_c:.2f} °C",
        "Evaporator minimum approach": f"{result.evaporator.minimum_approach_k:.2f} K",
        "Condenser minimum approach": f"{result.condenser.minimum_approach_k:.2f} K",
        "Energy balance error (Qc − Qe − W)": f"{c.energy_balance_error_w:.2e} W",
    }
    st.dataframe(pd.DataFrame({"Quantity": rows.keys(), "Value": rows.values()}), hide_index=True, **STRETCH)
    with st.expander("Cycle on the pressure-enthalpy diagram", expanded=True):
        show_figure(cycle_ph_figure(c))
    with st.expander("Heat-exchanger temperature profiles"):
        show_figure(hx_profiles_figure(result))
    with st.expander("State points"):
        st.dataframe(pd.DataFrame([s.as_display_dict() for s in c.states.values()]).round(3),
                     hide_index=True, **STRETCH)
    if not map_driven:
        st.download_button("Download traceable case result (JSON)", pd.Series(result.traceable_summary()).to_json(indent=2),
                           "coupled_case_result.json", "application/json")


solver_tab, sweep_tab, evaporator_tab, coil_tab, map_tab, validation_tab, fidelity_tab = st.tabs(
    ["1 · Coupled Solver", "2 · Design Sweep", "3 · Evaporator", "4 · Coil Geometry", "5 · Compressor Map", "6 · Test Validation", "7 · Model Fidelity"])

# ------------------------------------------------------------ coupled solver
with solver_tab:
    inp, result_panel = st.columns([0.36, 0.64])
    with inp:
        st.subheader("Operating conditions")
        fluid = st.selectbox("Refrigerant", REFRIGERANTS,
                             help="Must match the refrigerant of any imported compressor map. R407C is a zeotropic blend; saturation temperatures are dew points.")
        source_in = st.number_input("Source inlet (°C)", value=12.0, min_value=-40.0, max_value=60.0,
                                    help="Temperature of the fluid being cooled (e.g. chilled-water return) entering the evaporator.")
        sink_in = st.number_input("Sink inlet (°C)", value=30.0, min_value=-20.0, max_value=70.0,
                                  help="Temperature of the cooling medium (e.g. condenser water or ambient air) entering the condenser.")
        mdot = st.number_input("Refrigerant mass flow (kg/s)", min_value=0.001, max_value=2.0, value=0.03, step=0.001, format="%.4f",
                               help="Set by the compressor. Ignored when a compressor map closes the operating point.")
        eta = st.slider("Compressor isentropic efficiency", .40, .90, .72, .01,
                        help="Ideal (isentropic) compression work divided by actual work.")
        source_c = st.number_input("Source capacity rate (W/K)", min_value=50.0, value=1000.0,
                                   help="Mass flow × specific heat of the source stream. 1000 W/K ≈ 0.24 kg/s of water.")
        sink_c = st.number_input("Sink capacity rate (W/K)", min_value=50.0, value=1250.0,
                                 help="Mass flow × specific heat of the sink stream.")
        evap_ua = st.number_input("Evaporator UA (W/K)", min_value=10.0, value=950.0,
                                  help="Overall heat-transfer coefficient × area. Estimate it in tab 4 from geometry.")
        cond_ua = st.number_input("Condenser UA (W/K)", min_value=10.0, value=1200.0)
        aux_power = st.number_input("Declared fan + pump + controls power (W)", min_value=0.0, value=180.0, step=10.0,
                                    help="Included explicitly in system COP. Keep the electrical boundary consistent with any compressor map.")
        with st.expander("Superheat, subcooling and map rating conditions"):
            superheat = st.number_input("Evaporator outlet superheat (K)", min_value=0.0, max_value=30.0, value=5.0, step=0.5)
            subcooling = st.number_input("Condenser outlet subcooling (K)", min_value=0.0, max_value=30.0, value=3.0, step=0.5)
            info = st.session_state.get("map_info", {})
            rating_sh = st.number_input("Map rating superheat (K)", min_value=0.0, max_value=30.0,
                                        value=float(info.get("rating_superheat_k", 5.0)), step=0.1,
                                        help="Superheat at which the compressor map was published. Map capacity is converted to mass flow at these rating conditions, then corrected to your operating superheat/subcooling.")
            rating_sc = st.number_input("Map rating subcooling (K)", min_value=0.0, max_value=30.0,
                                        value=float(info.get("rating_subcooling_k", 3.0)), step=0.1)
        case_name = st.text_input("Case name / revision", value="5 kW water chiller - preliminary case")
        use_map = st.checkbox("Close operating point with imported compressor map", disabled="compressor_fit" not in st.session_state,
                              help="Fit a compressor map in tab 5 first. The operating point must stay inside the map envelope.")
        if use_map and info.get("refrigerant") and info["refrigerant"] != fluid:
            st.warning(f"The fitted map is for {info['refrigerant']}, but {fluid} is selected.")
        run = st.button("Run coupled CoolProp simulation", type="primary", **STRETCH)

    current = dict(fluid=fluid, source_in=source_in, sink_in=sink_in, mdot=mdot, eta=eta, source_c=source_c, sink_c=sink_c,
                   evap_ua=evap_ua, cond_ua=cond_ua, aux=aux_power, sh=superheat, sc=subcooling, rsh=rating_sh, rsc=rating_sc,
                   use_map=use_map, speed=st.session_state.get("map_speed", 3600.0))
    with result_panel:
        if run:
            st.session_state.pop("solver_error", None)
            try:
                with st.spinner("Solving the coupled refrigerant cycle..."):
                    if use_map:
                        st.session_state["solver_result"] = (True, solve_map_coupled_system(MapCoupledInputs(
                            fluid=fluid, speed_rpm=current["speed"], superheat_k=superheat, subcooling_k=subcooling,
                            map_rating_superheat_k=rating_sh, map_rating_subcooling_k=rating_sc,
                            source_inlet_c=source_in, source_capacity_rate_w_k=source_c, evaporator_ua_w_k=evap_ua,
                            sink_inlet_c=sink_in, sink_capacity_rate_w_k=sink_c, condenser_ua_w_k=cond_ua,
                            auxiliary_power_w=aux_power), st.session_state["compressor_fit"]))
                    else:
                        st.session_state["solver_result"] = (False, solve_coupled_system(CoupledInputs(
                            fluid=fluid, refrigerant_mass_flow_kg_s=mdot, compressor_isentropic_efficiency=eta,
                            superheat_k=superheat, subcooling_k=subcooling,
                            source_inlet_c=source_in, source_capacity_rate_w_k=source_c, evaporator_ua_w_k=evap_ua,
                            sink_inlet_c=sink_in, sink_capacity_rate_w_k=sink_c, condenser_ua_w_k=cond_ua,
                            auxiliary_power_w=aux_power, case_name=case_name)))
                    st.session_state["solver_inputs"] = current
            except Exception as exc:  # noqa: BLE001 - surface any model error to the user
                st.session_state.pop("solver_result", None)
                st.session_state["solver_inputs"] = current
                st.session_state["solver_error"] = str(exc)
        if err := st.session_state.get("solver_error"):
            stale_note("solver_inputs", current)
            st.error(f"No valid operating point: {err}")
        elif stored := st.session_state.get("solver_result"):
            stale_note("solver_inputs", current)
            show_operating_point(stored[1], map_driven=stored[0])
        else:
            st.info("Set the operating conditions and click **Run coupled CoolProp simulation**. The solver searches for the "
                    "evaporating and condensing temperatures at which the UA each heat exchanger needs equals the UA it has, "
                    "then reports capacity, COP, and the cycle on a pressure-enthalpy diagram.")

# -------------------------------------------------------------- design sweep
with sweep_tab:
    st.subheader("Design sweep - how does the system respond?")
    st.caption("Uses the operating conditions from tab 1 as the base case and re-solves the full coupled cycle at every point, "
               "so the evaporating and condensing temperatures move with the design change.")
    sweep_defaults = {"evaporator_ua_w_k": (300.0, 3000.0), "condenser_ua_w_k": (400.0, 4000.0), "sink_inlet_c": (20.0, 45.0),
                      "source_inlet_c": (6.0, 20.0), "refrigerant_mass_flow_kg_s": (0.01, 0.05),
                      "compressor_isentropic_efficiency": (0.5, 0.85), "sink_capacity_rate_w_k": (400.0, 4000.0),
                      "source_capacity_rate_w_k": (400.0, 4000.0)}
    left, right = st.columns([.32, .68])
    with left:
        field = st.selectbox("Variable to sweep", list(SWEEPABLE), format_func=lambda k: f"{SWEEPABLE[k][0]} ({SWEEPABLE[k][1]})")
        lo_d, hi_d = sweep_defaults[field]
        step = 0.001 if field == "refrigerant_mass_flow_kg_s" else (0.01 if field == "compressor_isentropic_efficiency" else 1.0)
        fmt = "%.3f" if step < 0.01 else ("%.2f" if step < 1 else "%.1f")
        lo = st.number_input("From", value=lo_d, step=step, format=fmt, key=f"sw_lo_{field}")
        hi = st.number_input("To", value=hi_d, step=step, format=fmt, key=f"sw_hi_{field}")
        n_points = st.slider("Points", 4, 15, 10)
        metric = st.selectbox("Headline metric", ["cycle_cop", "system_cop", "cooling_capacity_kw", "compressor_power_kw"],
                              format_func=lambda m: {"cycle_cop": "Cycle COP", "system_cop": "System COP",
                                                     "cooling_capacity_kw": "Cooling capacity (kW)", "compressor_power_kw": "Compressor power (kW)"}[m])
        go = st.button("Run sweep", type="primary")
    with right:
        if go:
            if hi <= lo:
                st.session_state["sweep"] = "'To' must be greater than 'From'."
            else:
                base = CoupledInputs(fluid=fluid, refrigerant_mass_flow_kg_s=mdot, compressor_isentropic_efficiency=eta,
                                     superheat_k=superheat, subcooling_k=subcooling, source_inlet_c=source_in,
                                     source_capacity_rate_w_k=source_c, evaporator_ua_w_k=evap_ua, sink_inlet_c=sink_in,
                                     sink_capacity_rate_w_k=sink_c, condenser_ua_w_k=cond_ua, auxiliary_power_w=aux_power)
                bar = st.progress(0.0, text="Solving sweep points...")
                df = sweep_coupled(base, field, np.linspace(lo, hi, n_points),
                                   progress=lambda i, n: bar.progress(i / n, text=f"Solved {i} of {n} points"))
                bar.empty()
                st.session_state["sweep"] = (field, metric, df)
        stored = st.session_state.get("sweep")
        if isinstance(stored, str):
            st.error(stored)
        elif stored:
            s_field, s_metric, df = stored
            label, unit = SWEEPABLE[s_field]
            ok = df[df["feasible"]]
            if ok.empty:
                st.error("No feasible operating point anywhere in this range.")
            else:
                names = {"cycle_cop": "Cycle COP", "system_cop": "System COP", "cooling_capacity_kw": "Cooling capacity (kW)",
                         "compressor_power_kw": "Compressor power (kW)"}
                fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
                axes[0].plot(ok[s_field], ok[s_metric], "o-", lw=2, color="#087f8c")
                axes[0].set_ylabel(names[s_metric])
                axes[0].set_title(f"{names[s_metric]} vs {label.lower()}")
                axes[1].plot(ok[s_field], ok["evaporating_c"], "o-", lw=2, label="evaporating")
                axes[1].plot(ok[s_field], ok["condensing_c"], "s-", lw=2, label="condensing")
                axes[1].set_ylabel("Saturation temperature (°C)")
                axes[1].set_title("Where the cycle settles")
                axes[1].legend()
                for ax in axes:
                    ax.set_xlabel(f"{label} ({unit})")
                    ax.grid(True, alpha=0.25)
                fig.tight_layout()
                show_figure(fig)
                if note := diminishing_returns(df, s_field, s_metric):
                    st.info(note)
                bad = df[~df["feasible"]]
                if not bad.empty:
                    st.warning(f"{len(bad)} point(s) have no valid operating point (edge of the envelope): "
                               + "; ".join(f"{v:.4g} {unit}" for v in bad[s_field]))
                st.dataframe(df.round(4), hide_index=True, **STRETCH)
                st.download_button("Download sweep (CSV)", df.to_csv(index=False).encode(), f"sweep_{s_field}.csv", "text/csv")
        else:
            st.info("Pick a variable and range, then run the sweep. Try **Evaporator UA**: capacity barely changes at fixed "
                    "refrigerant flow, but COP rises because the evaporating temperature climbs toward the source temperature, "
                    "with clear diminishing returns.")

# --------------------------------------------------------------- evaporator
with evaporator_tab:
    st.subheader("Evaporator task, rating, and refrigerant state checks")
    st.caption("Required cooling is not an equipment rating; constant-saturation UA is a bounded screening model; the refrigerant panel is a property/energy check, not a geometry prediction.")
    task_tab, screen_tab, states_tab = st.tabs(["Cooling requirement", "UA screening", "Refrigerant states"])
    with task_tab:
        left, right = st.columns([.42, .58])
        with left:
            secondary_kind = st.selectbox("Secondary-fluid type", ["Water", "Water-glycol (MEG 30%)", "Air (dry-air screening)", "Brine / process fluid", "Custom fluid"],
                                          help="The fluid on the non-refrigerant side of the evaporator.")
            coolprop_fluid = {"Water": "Water", "Water-glycol (MEG 30%)": "INCOMP::MEG-30%", "Air (dry-air screening)": "Air"}.get(secondary_kind)
            task_inlet = st.number_input("Secondary-fluid inlet temperature (°C)", value=12.0 if "Air" not in secondary_kind else 27.0, key=f"task_inlet_{secondary_kind}")
            task_target = st.number_input("Desired secondary-fluid outlet temperature (°C)", value=7.0 if "Air" not in secondary_kind else 15.0, key=f"task_target_{secondary_kind}")
            task_mdot = st.number_input("Secondary-fluid mass flow (kg/s)", min_value=.0001, value=(.4 if "Air" in secondary_kind else 1000/4180),
                                        format="%.6f", key=f"task_mdot_{secondary_kind}")
            if coolprop_fluid:
                use_cp = st.checkbox("Use constant specific heat instead of CoolProp (explicit approximation)", value=False, key=f"use_cp_{secondary_kind}")
            else:
                use_cp = True
                st.caption("CoolProp has no model for a generic brine/process fluid, so a measured or documented specific heat is required.")
            default_cp = 1006.0 if "Air" in secondary_kind else (3750.0 if "glycol" in secondary_kind else 4180.0)
            cp = st.number_input("Specific heat (J/kg-K)", min_value=1.0, value=default_cp, disabled=not use_cp, key=f"cp_{secondary_kind}")
            task_inputs = dict(kind=secondary_kind, tin=task_inlet, tout=task_target, m=task_mdot, use_cp=use_cp, cp=cp)
            if st.button("Calculate required cooling", type="primary"):
                try:
                    st.session_state["task_result"] = (liquid_cooling_task(fluid=coolprop_fluid or "Water", mass_flow_kg_s=task_mdot,
                        inlet_c=task_inlet, target_outlet_c=task_target, cp_j_kg_k=cp if use_cp else None), secondary_kind)
                    st.session_state["task_inputs"] = task_inputs
                except Exception as exc:  # noqa: BLE001
                    st.session_state["task_result"] = str(exc)
        with right:
            stored = st.session_state.get("task_result")
            if isinstance(stored, str):
                st.error(stored)
            elif stored:
                task, kind = stored
                stale_note("task_inputs", task_inputs)
                st.metric("Required cooling duty", f"{task.required_duty_w/1000:.3f} kW")
                st.write(f"Method: **{task.method}** ({kind}). This is the cooling the secondary fluid needs; it does not prove that a selected evaporator or compressor can supply it.")
                if "Air" in kind:
                    st.warning("Dry-air sensible screening only. If the coil is below the air dew point, use a wet-coil model with humidity and condensate before using this for design.")
                elif kind in {"Brine / process fluid", "Custom fluid"}:
                    st.warning("Uses your specific heat. Confirm it from a datasheet at the operating temperature and concentration.")
    with screen_tab:
        left, right = st.columns([.42, .58])
        with left:
            scr_ua = st.number_input("Evaporator UA (W/K)", min_value=1.0, value=950.0, key="screen_ua")
            scr_c = st.number_input("Secondary-fluid capacity rate (W/K)", min_value=1.0, value=1000.0, key="screen_c",
                                    help="Mass flow × specific heat of the fluid being cooled.")
            scr_inlet = st.number_input("Secondary-fluid inlet (°C)", value=12.0, key="screen_inlet")
            scr_sat = st.number_input("Assumed refrigerant saturation temperature (°C)", value=4.0, key="screen_sat")
            screen_inputs = dict(ua=scr_ua, c=scr_c, tin=scr_inlet, tsat=scr_sat)
            if st.button("Run UA screening", type="primary"):
                try:
                    st.session_state["screen_result"] = constant_saturation_evaporator(ua_w_k=scr_ua, capacity_rate_w_k=scr_c,
                        secondary_inlet_c=scr_inlet, saturation_temperature_c=scr_sat)
                    st.session_state["screen_inputs"] = screen_inputs
                except Exception as exc:  # noqa: BLE001
                    st.session_state["screen_result"] = str(exc)
        with right:
            screen = st.session_state.get("screen_result")
            if isinstance(screen, str):
                st.error(screen)
            elif screen:
                stale_note("screen_inputs", screen_inputs)
                a, b, c = st.columns(3)
                a.metric("Delivered duty", f"{screen.duty_w/1000:.3f} kW")
                b.metric("Secondary-fluid outlet", f"{screen.secondary_outlet_c:.3f} °C")
                c.metric("Effectiveness", f"{screen.effectiveness:.4f}", help="ε = 1 − exp(−NTU), NTU = UA / C, because the boiling refrigerant stays at constant temperature.")
                st.dataframe(pd.DataFrame({"Check": ["NTU", "LMTD", "UA back-calculated from Q/LMTD", "Assumptions"],
                    "Value": [f"{screen.ua_w_k/screen.capacity_rate_w_k:.3f}", f"{screen.lmtd_k:.3f} K", f"{screen.ua_check_w_k:.2f} W/K",
                              "Fixed saturation temperature; no refrigerant pressure loss or superheat zone"]}),
                    hide_index=True, **STRETCH)
    with states_tab:
        left, right = st.columns([.42, .58])
        with left:
            ref_fluid = st.selectbox("Pure refrigerant", ["R134a", "R290", "R32", "R1234yf"], key="evap_ref_fluid",
                                     help="Blends such as R407C need a dew/bubble-aware model and are rejected here.")
            ref_sat = st.number_input("Saturation temperature (°C)", value=5.0, key="ref_sat")
            ref_x = st.number_input("Refrigerant inlet quality (-)", min_value=0.0, max_value=1.0, value=.2, step=.05,
                                    help="Vapor mass fraction leaving the expansion valve, typically 0.2-0.35.")
            ref_out = st.number_input("Refrigerant outlet temperature (°C)", value=10.0, key="ref_out",
                                      help="Must be at or above the saturation temperature (no pressure drop is modelled).")
            ref_mdot = st.number_input("Refrigerant mass flow (kg/s)", min_value=.00001, value=.004, format="%.5f", key="ref_mdot")
            states_inputs = dict(f=ref_fluid, tsat=ref_sat, x=ref_x, tout=ref_out, m=ref_mdot)
            if st.button("Check refrigerant energy balance", type="primary"):
                try:
                    st.session_state["states_result"] = refrigerant_evaporator_states(fluid=ref_fluid, saturation_temperature_c=ref_sat,
                        inlet_quality=ref_x, outlet_temperature_c=ref_out, refrigerant_mass_flow_kg_s=ref_mdot)
                    st.session_state["states_inputs"] = states_inputs
                except Exception as exc:  # noqa: BLE001
                    st.session_state["states_result"] = str(exc)
        with right:
            props = st.session_state.get("states_result")
            if isinstance(props, str):
                st.error(props)
            elif props:
                stale_note("states_inputs", states_inputs)
                a, b, c = st.columns(3)
                a.metric("Saturation pressure", f"{props.saturation_pressure_pa/1000:.2f} kPa(a)")
                b.metric("Refrigerant enthalpy duty", f"{props.duty_w:.1f} W")
                c.metric("Outlet superheat", f"{props.superheat_k:.2f} K")
                st.write(f"Outlet status: **{props.outlet_phase}**. Inlet enthalpy {props.inlet_enthalpy_j_kg/1000:.1f} kJ/kg, "
                         f"outlet {props.outlet_enthalpy_j_kg/1000:.1f} kJ/kg.")
    st.warning("Current physical limit: no validated plate/finned-tube boiling correlation with refrigerant pressure-drop closure yet. Do not treat UA screening or property checks as a geometry rating.")

# ------------------------------------------------------------ coil geometry
with coil_tab:
    st.subheader("Tube-bundle evaporator or condenser")
    st.caption("Enter drawing or datasheet values. The refrigerant-side h must be measured or come from a documented correlation valid for this refrigerant, geometry, and flow regime.")
    a, b, c, d = st.columns(4)
    with a: side = st.selectbox("Component", ["Evaporator", "Condenser"])
    with b: di = st.number_input("Tube inner diameter (mm)", 1.0, 100.0, value=8.0)
    with c: do = st.number_input("Tube outer diameter (mm)", 1.1, 110.0, value=10.0)
    with d: length = st.number_input("Length per tube (m)", .1, 50.0, value=4.0)
    e, f, g, h = st.columns(4)
    with e: count = st.number_input("Tube count", 1, 10000, value=24)
    with f: h_ref = st.number_input("Refrigerant-side h (W/m²·K)", 1.0, value=1600.0, help="Inside the tubes.")
    with g: h_sec = st.number_input("Secondary-side h (W/m²·K)", 1.0, value=120.0, help="Outside the tubes. ~50-150 for air, ~1000-5000 for water.")
    with h: fouling = st.number_input("Outer fouling resistance (m²·K/W)", 0.0, 0.01, value=0.0, format="%.6f")
    status = st.selectbox("Correlation / coefficient provenance", ["Measured or calibrated", "Documented single-phase correlation", "Documented boiling correlation", "Documented condensation correlation", "Placeholder — not suitable for validation"])
    try:
        ua = estimate_tube_bundle_ua(TubeBundleGeometry(di/1000, do/1000, length, int(count), fouling_resistance_outer_m2k_w=fouling), h_ref, h_sec)
        i, j, k = st.columns(3)
        i.metric(f"Estimated {side} UA", f"{ua.ua_w_k:.0f} W/K")
        j.metric("Tube areas (inner / outer)", f"{ua.inner_area_m2:.2f} / {ua.outer_area_m2:.2f} m²")
        k.metric("Overall U, outer basis", f"{ua.overall_u_outer_w_m2_k:.1f} W/m²·K")
        resistances = {"Inner film": ua.resistance_inner_k_w, "Wall": ua.resistance_wall_k_w, "Outer film": ua.resistance_outer_k_w, "Fouling": ua.resistance_fouling_k_w}
        total = sum(resistances.values())
        st.dataframe(pd.DataFrame({"Resistance": resistances.keys(), "K/W": [f"{v:.3e}" for v in resistances.values()],
                                   "Share of total": [f"{v/total:.1%}" for v in resistances.values()]}), hide_index=True, **STRETCH)
        st.caption(f"Coefficient provenance: {status}. The largest resistance controls UA - improve that side first.")
    except ValueError as exc:
        st.error(str(exc))
    st.divider()
    st.subheader("Plate evaporator geometry path - preliminary")
    st.caption("Counterflow plate geometry with secondary-side heat transfer (smooth parallel-plate channel: Nu = 7.54 laminar, Gnielinski turbulent, interpolated in transition) and pressure drop (f = 96/Re laminar). The refrigerant boiling coefficient is prescribed.")
    p1, p2, p3, p4 = st.columns(4)
    with p1: plate_length = st.number_input("Plate active length (m)", .05, 3.0, value=.50)
    with p2: plate_width = st.number_input("Plate active width (m)", .03, 1.5, value=.12)
    with p3: plates = st.number_input("Plate count", 4, 500, value=30)
    with p4: gap_mm = st.number_input("Channel gap (mm)", .5, 20.0, value=2.5)
    p5, p6, p7, p8 = st.columns(4)
    with p5: enlargement = st.number_input("Area enlargement factor", 1.0, 2.0, value=1.15)
    with p6: plate_flow = st.number_input("Secondary-fluid flow (kg/s)", .001, value=.24)
    with p7: plate_tin = st.number_input("Secondary-fluid inlet (°C)", value=12.0, key="plate_tin")
    with p8: plate_tsat = st.number_input("Refrigerant saturation (°C)", value=4.0, key="plate_tsat")
    q1, q2 = st.columns(2)
    with q1: plate_href = st.number_input("Prescribed refrigerant boiling h (W/m²-K)", 1.0, value=1800.0, help="Enter only a measured/calibrated or documented correlation result.")
    with q2: plate_fluid_label = st.selectbox("Secondary fluid", ["Water", "Water-glycol (MEG 30%)"], key="plate_fluid")
    plate_fluid = "Water" if plate_fluid_label == "Water" else "INCOMP::MEG-30%"
    plate_inputs = dict(L=plate_length, W=plate_width, n=plates, gap=gap_mm, phi=enlargement, m=plate_flow, tin=plate_tin, tsat=plate_tsat, h=plate_href, f=plate_fluid)
    if st.button("Run plate-evaporator screening", type="primary"):
        try:
            st.session_state["plate_result"] = plate_evaporator_screening(
                geometry=PlateEvaporatorGeometry(plate_length, plate_width, int(plates), gap_mm/1000, enlargement_factor=enlargement),
                secondary_fluid=plate_fluid, secondary_mass_flow_kg_s=plate_flow, secondary_inlet_c=plate_tin,
                refrigerant_saturation_c=plate_tsat, refrigerant_h_w_m2_k=plate_href)
            st.session_state["plate_inputs"] = plate_inputs
        except Exception as exc:  # noqa: BLE001
            st.session_state["plate_result"] = str(exc)
    plate = st.session_state.get("plate_result")
    if isinstance(plate, str):
        st.error(plate)
    elif plate:
        stale_note("plate_inputs", plate_inputs)
        x, y, z, w = st.columns(4)
        x.metric("Geometry UA", f"{plate.ua_w_k:.0f} W/K")
        y.metric("Delivered duty", f"{plate.screening.duty_w/1000:.2f} kW")
        z.metric("Secondary outlet", f"{plate.screening.secondary_outlet_c:.2f} °C")
        w.metric("Secondary-side Δp", f"{plate.secondary_pressure_drop_pa/1000:.2f} kPa")
        st.info(f"Secondary side: Re = {plate.secondary_reynolds:.0f} ({plate.secondary_regime}); Nu = {plate.secondary_nusselt:.2f}; "
                f"h = {plate.secondary_h_w_m2_k:.0f} W/m²-K. {plate.correlation_status}")

# ----------------------------------------------------------- compressor map
with map_tab:
    st.subheader("Compressor-map import, calibration, and held-out validation")
    st.write("Export the manufacturer map or controlled compressor-test data to CSV. Capacity and power must both be in **W**; "
             "rows that would beat the Carnot limit are rejected (usually a Btu/h vs W mix-up).")
    c1, c2 = st.columns(2)
    with c1:
        download_csv(COMPRESSOR_MAP_COLUMNS, [{"evaporating_temperature_c": -10, "condensing_temperature_c": 35, "speed_rpm": 3600, "capacity_w": 3900, "power_w": 1100},
                                              {"evaporating_temperature_c": -5, "condensing_temperature_c": 40, "speed_rpm": 3600, "capacity_w": 4500, "power_w": 1280}],
                     "Download compressor-map template")
    with c2:
        if EXAMPLE_MAP.exists() and st.button("Load bundled Copeland R407C example"):
            st.session_state["map_df"] = pd.read_csv(EXAMPLE_MAP)
            st.session_state["map_info"] = EXAMPLE_MAP_INFO
    uploaded = st.file_uploader("Upload compressor-map CSV", type="csv")
    if uploaded is not None and st.session_state.get("map_upload_id") != uploaded.file_id:
        st.session_state["map_upload_id"] = uploaded.file_id
        try:
            st.session_state["map_df"] = pd.read_csv(uploaded)
            st.session_state["map_info"] = {"name": uploaded.name}
        except Exception as exc:  # noqa: BLE001
            st.error(f"Could not read CSV: {exc}")
    st.number_input("Compressor speed for map closure (rpm)", 100.0, value=3600.0, step=100.0, key="map_speed")
    fraction = st.slider("Calibration share", .50, .90, .75, .05, help="Fraction of rows used to fit; the rest are held out to measure prediction error.")
    seed = st.number_input("Reproducible split seed", 1, value=2026, step=1)
    raw = st.session_state.get("map_df")
    if raw is None:
        st.info("Required CSV columns: " + ", ".join(COMPRESSOR_MAP_COLUMNS) + ". A minimum of 14 distinct map points is required.")
    else:
        info = st.session_state.get("map_info", {})
        try:
            data = validate_compressor_map(raw)
            digest = hashlib.sha1(pd.util.hash_pandas_object(data, index=False).values.tobytes()).hexdigest()
            if st.session_state.get("fit_digest") not in (None, digest):
                st.session_state.pop("compressor_fit", None)  # never keep a fit from a different map
            st.success(f"Validated {len(data)} distinct map points from {info.get('name', 'uploaded file')}"
                       + (f" - refrigerant {info['refrigerant']}." if info.get("refrigerant") else "."))
            st.dataframe(data, **STRETCH, height=240)
            if st.button("Fit map and score held-out rows", type="primary"):
                st.session_state["compressor_fit"] = fit_compressor_map(data, fraction, int(seed))
                st.session_state["fit_digest"] = digest
            if fit := st.session_state.get("compressor_fit"):
                st.markdown("**Held-out prediction error** (rows the fit never saw)")
                x, y = st.columns(2)
                with x:
                    st.caption("Capacity")
                    st.dataframe(metric_table(fit.capacity_validation), hide_index=True, **STRETCH)
                with y:
                    st.caption("Electrical power")
                    st.dataframe(metric_table(fit.power_validation), hide_index=True, **STRETCH)
                st.warning(f"Fit uses {fit.calibration_rows} calibration and {fit.validation_rows} held-out rows. Envelope: evaporating "
                           f"{fit.evaporating_range_c[0]:.1f} to {fit.evaporating_range_c[1]:.1f} °C, condensing {fit.condensing_range_c[0]:.1f} to "
                           f"{fit.condensing_range_c[1]:.1f} °C, {fit.speed_range_rpm[0]:.0f}-{fit.speed_range_rpm[1]:.0f} rpm. Extrapolation is disabled.")
                if info.get("rating_superheat_k") is not None:
                    st.caption(f"Map rating conditions: {info['rating_superheat_k']} K superheat, {info['rating_subcooling_k']} K subcooling "
                               "(pre-filled in the Coupled Solver).")
        except Exception as exc:  # noqa: BLE001
            st.session_state.pop("compressor_fit", None)
            st.error(f"Map rejected: {exc}")

# ---------------------------------------------------------- test validation
with validation_tab:
    st.subheader("Manufacturer datasheet check")
    st.write("The Copeland ZR34K3E-TFD datasheet publishes refrigerant **mass flow** and **isentropic efficiency** at 57 rating points. "
             "Neither is an input to this model, so they independently test its thermodynamics: the model predicts mass flow from "
             "published capacity (m = C / (h1 − h4) at the sheet's 11.1 K superheat / 8.3 K subcooling, dew-point convention) and "
             "isentropic efficiency from published mass flow and power (η = M (h2s − h1) / P).")
    if EXAMPLE_DATASHEET.exists() and st.button("Run datasheet check", type="primary"):
        with st.spinner("Evaluating 57 rating points with CoolProp..."):
            st.session_state["datasheet_check"] = validate_against_datasheet(
                pd.read_csv(EXAMPLE_DATASHEET), fluid="R407C", rating_superheat_k=20 / 1.8, rating_subcooling_k=15 / 1.8)
    if check := st.session_state.get("datasheet_check"):
        t = check.table
        a, b, c, d = st.columns(4)
        a.metric("Mass flow: mean abs. error", f"{check.mass_flow.mape_percent:.2f} %")
        b.metric("Mass flow: worst point", f"{check.mass_flow.maximum_absolute_percent_error:.2f} %")
        c.metric("Isentropic eff.: mean abs. error", f"{t.isentropic_eff_error_points.abs().mean():.2f} pts")
        d.metric("Isentropic eff.: worst point", f"{t.isentropic_eff_error_points.abs().max():.2f} pts")
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
        for ax, (x, y, title, unit) in zip(axes, [
                (t.mass_flow_sheet_kg_s * 1000, t.mass_flow_model_kg_s * 1000, "Refrigerant mass flow", "g/s"),
                (t.isentropic_eff_sheet_pct, t.isentropic_eff_model_pct, "Isentropic efficiency", "%")]):
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="perfect agreement")
            ax.scatter(x, y, c=t.condensing_c, cmap="viridis", s=28, zorder=3)
            ax.set_xlabel(f"Copeland datasheet ({unit})")
            ax.set_ylabel(f"This model ({unit})")
            ax.set_title(title)
            ax.grid(True, alpha=0.25)
            ax.legend(loc="upper left")
        fig.tight_layout()
        show_figure(fig)
        st.caption("Points coloured by condensing temperature. Published values are nominal (±5%) and rounded, so sub-1% "
                   "agreement is at the resolution of the sheet. This checks properties and state-point conventions, "
                   "not the heat-exchanger UA models.")
        with st.expander("Point-by-point table"):
            st.dataframe(t.round(4), hide_index=True, **STRETCH)
    st.divider()
    st.subheader("Physical test-point validation")
    st.write("Upload predictions generated from a frozen model configuration and measurements from independent controlled test points. Do not score the same points used to calibrate the model.")
    download_csv(TEST_POINT_COLUMNS, [{"test_id": "TP-01", "predicted_capacity_w": 4500, "measured_capacity_w": 4380, "predicted_power_w": 1280, "measured_power_w": 1310}], "Download test-comparison template")
    uploaded_tests = st.file_uploader("Upload prediction-versus-measurement CSV", type="csv", key="tests")
    if uploaded_tests:
        try:
            table = pd.read_csv(uploaded_tests)
            scores = score_test_points(table)
            st.dataframe(table, **STRETCH)
            a, b = st.columns(2)
            with a:
                st.caption("Capacity")
                st.dataframe(metric_table(scores["capacity"]), hide_index=True, **STRETCH)
            with b:
                st.caption("Electrical power")
                st.dataframe(metric_table(scores["power"]), hide_index=True, **STRETCH)
            st.warning("Set acceptance limits before examining results. Retain raw data, instrument calibration, uncertainty budget, conditions, software commit, and reviewer approval.")
        except Exception as exc:  # noqa: BLE001
            st.error(f"Test data rejected: {exc}")
    st.markdown("**Evidence to retain:** product configuration, manufacturer-map source, refrigerant, geometry drawings, correlation reference and validity range, charge, controls, calibrated instruments, raw data, and independent review.")

# ----------------------------------------------------------- model fidelity
with fidelity_tab:
    st.subheader("Model fidelity and evidence status")
    st.caption("This page separates what is calculated, what is calibrated, and what still needs independent evidence. It does not invent a confidence percentage.")
    st.markdown("""
| Level | What the model can answer | Status |
|---|---|---|
| 1 - Screening | Cooling requirement, UA-based duty, LMTD/effectiveness, property states | Implemented; unit-tested against hand calculations |
| 2 - Component prediction | Finite-UA coupled cycle, compressor-map closure, preliminary plate geometry | Implemented; validity depends on supplied UA, map, and correlations |
| 3 - Validated hardware model | Error against independent held-out measurements with uncertainty and configuration control | Workflow implemented; needs measured test data |
""")
    st.markdown("### Current evidence checklist")
    col_a, col_b = st.columns(2)
    with col_a:
        st.checkbox("Manufacturer compressor-map rows imported and envelope checked", value="compressor_fit" in st.session_state, disabled=True)
        st.checkbox("Independent hold-out map rows scored", value="compressor_fit" in st.session_state, disabled=True)
        st.checkbox("Selected evaporator geometry documented", value=False, disabled=True)
    with col_b:
        st.checkbox("Boiling and two-phase pressure-drop correlations selected for geometry", value=False, disabled=True)
        st.checkbox("Independent physical test points compared", value=False, disabled=True)
        st.checkbox("Uncertainty budget and configuration record retained", value=False, disabled=True)
    st.markdown("### Known limitations")
    st.markdown("""
- Steady state only; no start-up, cycling, or control dynamics.
- No refrigerant-side pressure drop; saturation pressure is constant through each exchanger.
- Plate secondary side uses smooth-channel correlations - real chevron plates transfer heat better, so this is conservative.
- Refrigerant boiling/condensing coefficients are user-supplied, not predicted.
- Air-side cooling is dry (sensible) only; no dehumidification.
""")
