"""Robust saturation-temperature closure shared by the coupled solvers.

Both coupled solvers look for the evaporating and condensing temperatures at
which the UA each exchanger *requires* equals the UA it *has*.  A fast bounded
least-squares solve is tried first.  If it stalls (for example on a flat
penalty plateau next to a temperature pinch), a bracketed nested root search
is used instead.  It relies on two physical monotonicities:

* raising the condensing temperature increases the condenser driving force,
  so the condenser's required UA falls;
* raising the evaporating temperature reduces the evaporator driving force,
  so the evaporator's required UA rises.

A pinch or temperature cross means "infinite UA required", which keeps both
functions monotone.  If no bracket exists the problem is physically
infeasible, and the error message says which exchanger is the limit.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import inf, isfinite, log
from typing import Callable

import numpy as np
from scipy.optimize import brentq, least_squares


class InfeasibleOperatingPoint(ValueError):
    """No steady operating point exists inside the allowed temperature bounds."""


@dataclass(frozen=True)
class ClosureResult:
    evaporating_c: float
    condensing_c: float
    method: str
    function_evaluations: int
    message: str


# required_ua(te_c, tc_c) -> (evaporator_required_ua, condenser_required_ua).
# Each value may be +inf for a pinch/cross or non-physical state.
RequiredUA = Callable[[float, float], tuple[float, float]]


def close_saturation_temperatures(
    required_ua: RequiredUA, *, evaporator_ua: float, condenser_ua: float,
    te_bounds: tuple[float, float], tc_bounds: tuple[float, float],
    initial: tuple[float, float], tolerance: float = 1e-5,
) -> ClosureResult:
    te_lo, te_hi = te_bounds
    tc_lo, tc_hi = tc_bounds
    if not (te_lo < te_hi and tc_lo < tc_hi):
        raise InfeasibleOperatingPoint(
            "The allowed temperature window is empty. Check that the source is colder than the sink, "
            "that the temperatures are inside the refrigerant's valid range, and (for map closure) "
            "that the compressor map envelope overlaps these conditions.")
    calls = 0
    cache: dict[tuple[float, float], tuple[float, float]] = {}

    def ua(te: float, tc: float) -> tuple[float, float]:
        nonlocal calls
        key = (round(te, 9), round(tc, 9))
        if key not in cache:
            calls += 1
            if tc <= te:
                cache[key] = (inf, inf)
            else:
                cache[key] = required_ua(te, tc)
        return cache[key]

    def rel(x: float, target: float) -> float:
        return (x - target) / target if isfinite(x) else 10.0

    # 1) Fast path: bounded least squares from a few sensible starting points.
    #    Cold, high-lift starts help when the warm default sits on a pinch plateau.
    starts = [initial,
              (te_hi - 0.35 * (te_hi - te_lo), tc_lo + 0.35 * (tc_hi - tc_lo)),
              (te_hi - 0.6 * (te_hi - te_lo), tc_lo + 0.6 * (tc_hi - tc_lo))]
    for start in starts:
        x0 = np.clip(np.array(start, float), [te_lo + 1e-6, tc_lo + 1e-6], [te_hi - 1e-6, tc_hi - 1e-6])
        try:
            sol = least_squares(
                lambda x: np.array([rel(ua(x[0], x[1])[0], evaporator_ua), rel(ua(x[0], x[1])[1], condenser_ua)]),
                x0, bounds=([te_lo, tc_lo], [te_hi, tc_hi]), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=120)
            ev, co = ua(float(sol.x[0]), float(sol.x[1]))
            if isfinite(ev) and isfinite(co) and max(abs(rel(ev, evaporator_ua)), abs(rel(co, condenser_ua))) < tolerance:
                return ClosureResult(float(sol.x[0]), float(sol.x[1]), "least_squares", calls, str(sol.message))
        except (ValueError, FloatingPointError):
            continue

    # 2) Fallback: nested bracketed root search on log(UA_required / UA_available).
    def log_ratio(x: float, target: float) -> float:
        return log(x / target) if isfinite(x) and x > 0 else inf

    def capped(value: float) -> float:
        # Required UA grows continuously toward a pinch, so a large finite cap
        # preserves the sign change that brackets the true root.
        return value if isfinite(value) else 50.0

    def condensing_for(te: float) -> float | None:
        """Condensing temperature that closes the condenser at this evaporating temperature."""
        f = lambda tc: capped(log_ratio(ua(te, tc)[1], condenser_ua))
        grid = np.linspace(tc_lo, tc_hi, 14)
        values = [f(float(tc)) for tc in grid]
        for i in range(len(grid) - 1):
            if values[i] > 0 >= values[i + 1]:
                return brentq(f, float(grid[i]), float(grid[i + 1]), xtol=1e-7)
        return None

    def evaporator_residual(te: float) -> float:
        tc = condensing_for(te)
        if tc is None:
            return float("nan")
        return log_ratio(ua(te, tc)[0], evaporator_ua)

    te_grid = np.linspace(te_hi, te_lo, 14)  # from warm (high UA needed) to cold
    g = [evaporator_residual(float(t)) for t in te_grid]
    any_condenser_closure = any(not np.isnan(v) for v in g)
    for i in range(len(te_grid) - 1):
        gi, gj = g[i], g[i + 1]
        if np.isnan(gi) or np.isnan(gj):
            continue
        if gi > 0 >= gj:
            a, b = float(te_grid[i + 1]), float(te_grid[i])  # a colder (g<=0), b warmer (g>0)

            def safe(te: float) -> float:
                v = evaporator_residual(te)
                if np.isnan(v):
                    raise InfeasibleOperatingPoint("Condenser closure lost during evaporator search.")
                return capped(v)

            try:
                te = brentq(safe, a, b, xtol=1e-7)
            except InfeasibleOperatingPoint:
                break
            tc = condensing_for(te)
            ev, co = ua(te, tc)
            if max(abs(rel(ev, evaporator_ua)), abs(rel(co, condenser_ua))) < max(tolerance, 1e-5):
                return ClosureResult(te, tc, "bracketed_fallback", calls, "Converged by bracketed nested root search.")

    if not any_condenser_closure:
        raise InfeasibleOperatingPoint(
            "The condenser cannot reject the heat within the allowed condensing-temperature range. "
            "Increase condenser UA or sink capacity rate, lower the sink inlet temperature, or reduce refrigerant flow.")
    raise InfeasibleOperatingPoint(
        "The evaporator cannot absorb the refrigerant-side duty without a temperature cross. "
        "Increase evaporator UA or source capacity rate, raise the source inlet temperature, or reduce refrigerant flow.")


def fluid_temperature_limits_c(fluid: str) -> tuple[float, float]:
    """(minimum usable temperature, critical temperature) in °C, with a clear error for unknown fluids."""
    from CoolProp.CoolProp import PropsSI

    try:
        t_crit = float(PropsSI("Tcrit", fluid)) - 273.15
        t_min = float(PropsSI("Tmin", fluid)) - 273.15
    except ValueError as exc:
        raise ValueError(f"CoolProp does not recognise refrigerant '{fluid}'.") from exc
    return t_min, t_crit


def cycle_temperature_bounds(fluid: str, source_inlet_c: float, sink_inlet_c: float) -> tuple[tuple[float, float], tuple[float, float]]:
    """Physically allowed evaporating/condensing temperature windows."""
    t_min, t_crit = fluid_temperature_limits_c(fluid)
    if source_inlet_c <= t_min + 5:
        raise ValueError(f"Source inlet {source_inlet_c:.1f} °C is too cold for {fluid}.")
    if sink_inlet_c >= t_crit - 5:
        raise ValueError(f"Sink inlet {sink_inlet_c:.1f} °C is too close to the {fluid} critical temperature ({t_crit:.1f} °C); "
                         "a subcritical condensing cycle is not possible.")
    te = (max(source_inlet_c - 40.0, t_min + 5.0), source_inlet_c - 0.2)
    tc = (sink_inlet_c + 0.2, min(sink_inlet_c + 50.0, t_crit - 3.0))
    return te, tc
