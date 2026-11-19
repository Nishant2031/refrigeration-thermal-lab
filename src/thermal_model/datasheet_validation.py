"""Independent check of the cycle thermodynamics against a manufacturer datasheet.

A compressor datasheet (AHRI 540 format) publishes, at each rating point,
capacity C, electrical power P, refrigerant mass flow M and isentropic
efficiency. Mass flow and efficiency are NOT inputs to this model, so they
make an independent test of its state-point construction:

* predicted mass flow   m = C / (h1 - h4), at the sheet's rating superheat and
  subcooling, with dew-point saturation temperatures;
* predicted efficiency  eta = M (h2s - h1) / P, using the datasheet mass flow
  and CoolProp's isentropic enthalpy rise.

Agreement verifies the property model, the dew/bubble convention, and the
superheat/subcooling handling. It is not a validation of the heat-exchanger
UA models, which need system test data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

from .cycle import CycleInputs, solve_cycle
from .data_workflow import ErrorMetrics, error_metrics

BTU_H_TO_W = 0.29307107
LB_H_TO_KG_S = 0.45359237 / 3600.0
DATASHEET_COLUMNS = ("evaporating_temperature_f", "condensing_temperature_f", "capacity_btu_h", "power_w",
                     "mass_flow_lb_h", "isentropic_efficiency_pct")


def f_to_c(t_f):
    return (np.asarray(t_f, float) - 32.0) / 1.8


@dataclass(frozen=True)
class DatasheetValidation:
    table: pd.DataFrame
    mass_flow: ErrorMetrics
    isentropic_efficiency: ErrorMetrics


def validate_against_datasheet(data: pd.DataFrame, *, fluid: str, rating_superheat_k: float,
                               rating_subcooling_k: float) -> DatasheetValidation:
    missing = [c for c in DATASHEET_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(f"Datasheet is missing columns: {', '.join(missing)}")
    rows = []
    for r in data.itertuples(index=False):
        te, tc = float(f_to_c(r.evaporating_temperature_f)), float(f_to_c(r.condensing_temperature_f))
        cyc = solve_cycle(CycleInputs(fluid=fluid, evaporating_temperature_c=te, condensing_temperature_c=tc,
                                      superheat_k=rating_superheat_k, subcooling_k=rating_subcooling_k,
                                      refrigerant_mass_flow_kg_s=1.0))
        h1, h4 = cyc.states["1"].enthalpy_j_kg, cyc.states["4"].enthalpy_j_kg
        h2s = cyc.state_2s.enthalpy_j_kg
        capacity_w = r.capacity_btu_h * BTU_H_TO_W
        m_sheet = r.mass_flow_lb_h * LB_H_TO_KG_S
        m_pred = capacity_w / (h1 - h4)
        eta_pred = m_sheet * (h2s - h1) / r.power_w * 100.0
        rows.append(dict(evaporating_c=te, condensing_c=tc, capacity_w=capacity_w, power_w=r.power_w,
                         mass_flow_sheet_kg_s=m_sheet, mass_flow_model_kg_s=m_pred,
                         mass_flow_error_pct=100 * (m_pred - m_sheet) / m_sheet,
                         isentropic_eff_sheet_pct=r.isentropic_efficiency_pct, isentropic_eff_model_pct=eta_pred,
                         isentropic_eff_error_points=eta_pred - r.isentropic_efficiency_pct))
    t = pd.DataFrame(rows)
    return DatasheetValidation(
        t, error_metrics(t.mass_flow_model_kg_s, t.mass_flow_sheet_kg_s),
        error_metrics(t.isentropic_eff_model_pct, t.isentropic_eff_sheet_pct))


def sheet_to_compressor_map(data: pd.DataFrame, speed_rpm: float = 3600.0) -> pd.DataFrame:
    """Convert a °F / Btu/h datasheet table into the app's SI compressor-map CSV schema."""
    return pd.DataFrame({
        "evaporating_temperature_c": f_to_c(data["evaporating_temperature_f"]).round(3),
        "condensing_temperature_c": f_to_c(data["condensing_temperature_f"]).round(3),
        "speed_rpm": speed_rpm,
        "capacity_w": (data["capacity_btu_h"] * BTU_H_TO_W).round(1),
        "power_w": data["power_w"],
        "capacity_btu_h": data["capacity_btu_h"],
    })
