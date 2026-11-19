"""One-parameter design sweeps of the coupled cycle.

A sweep re-solves the full coupled operating point at each value, so it shows
how the *system* responds - including the saturation temperatures moving -
rather than a single component in isolation.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Callable

import numpy as np
import pandas as pd

from .coupled import CoupledInputs, solve_coupled_system

SWEEPABLE = {
    "evaporator_ua_w_k": ("Evaporator UA", "W/K"),
    "condenser_ua_w_k": ("Condenser UA", "W/K"),
    "sink_inlet_c": ("Sink inlet temperature", "°C"),
    "source_inlet_c": ("Source inlet temperature", "°C"),
    "refrigerant_mass_flow_kg_s": ("Refrigerant mass flow", "kg/s"),
    "compressor_isentropic_efficiency": ("Compressor isentropic efficiency", "-"),
    "sink_capacity_rate_w_k": ("Sink capacity rate", "W/K"),
    "source_capacity_rate_w_k": ("Source capacity rate", "W/K"),
}


def sweep_coupled(base: CoupledInputs, field: str, values,
                  progress: Callable[[int, int], None] | None = None) -> pd.DataFrame:
    """Solve the coupled cycle for each value of ``field``.

    Infeasible points are kept as rows with ``feasible=False`` and the reason,
    so a sweep never silently drops the edge of the operating envelope.
    """
    if field not in SWEEPABLE:
        raise ValueError(f"Cannot sweep '{field}'. Choose one of: {', '.join(SWEEPABLE)}")
    values = [float(v) for v in values]
    if len(values) < 2:
        raise ValueError("A sweep needs at least two values.")
    rows = []
    for i, value in enumerate(values):
        row = {field: value}
        try:
            r = solve_coupled_system(replace(base, **{field: value}))
            c = r.cycle
            row.update(feasible=r.solver_success, reason="" if r.solver_success else "did not converge",
                       cooling_capacity_kw=c.evaporator_heat_w / 1000, compressor_power_kw=c.compressor_power_w / 1000,
                       cycle_cop=c.cop_cooling, system_cop=r.system_cop_cooling,
                       evaporating_c=c.inputs.evaporating_temperature_c, condensing_c=c.inputs.condensing_temperature_c,
                       second_law_efficiency=c.cop_cooling / c.carnot_cop_cooling)
        except ValueError as exc:
            row.update(feasible=False, reason=str(exc))
        rows.append(row)
        if progress:
            progress(i + 1, len(values))
    return pd.DataFrame(rows)


def diminishing_returns(df: pd.DataFrame, field: str, metric: str = "cooling_capacity_kw") -> str | None:
    """Plain-language comparison of the metric's slope over the first vs last third of the sweep."""
    ok = df[df["feasible"]].sort_values(field)
    if len(ok) < 4:
        return None
    x, y = ok[field].to_numpy(float), ok[metric].to_numpy(float)
    k = max(2, len(ok) // 3)
    first = np.polyfit(x[:k], y[:k], 1)[0]
    last = np.polyfit(x[-k:], y[-k:], 1)[0]
    label, unit = SWEEPABLE[field]
    m_label = metric.replace("_kw", "").replace("_", " ")
    if abs(first) < 1e-12:
        return None
    ratio = last / first
    if ratio < 0:
        trend = "reverses direction across the range"
    elif ratio < 0.7:
        trend = f"weakens to {ratio:.0%} of its starting strength - diminishing returns"
    elif ratio > 1.3:
        trend = f"strengthens to {ratio:.1f}x its starting strength"
    else:
        trend = "stays roughly constant"
    return (f"{m_label.capitalize()} changes by {first:.3g} per {unit} at the low end of {label} and "
            f"{last:.3g} per {unit} at the high end: the effect {trend}.")
