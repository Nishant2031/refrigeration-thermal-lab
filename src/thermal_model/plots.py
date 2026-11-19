"""Plotting functions kept separate from the numerical model."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import PropsSI

from .coupled import CoupledResult
from .cycle import CycleInputs, CycleResult, solve_cycle


def plot_cycle_ph(result: CycleResult, path: str | Path) -> None:
    """Pressure-enthalpy diagram with a saturation dome, saved to ``path``."""
    fig = cycle_ph_figure(result)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def cycle_ph_figure(result: CycleResult):
    """Return a matplotlib Figure of the cycle on a pressure-enthalpy diagram."""
    f = result.inputs.fluid
    tmin = PropsSI("Ttriple", f) + 1.0
    tcrit = PropsSI("Tcrit", f) - 0.2
    temps = np.linspace(tmin, tcrit, 300)
    hf = np.array([PropsSI("Hmass", "T", t, "Q", 0, f) for t in temps]) / 1e3
    hg = np.array([PropsSI("Hmass", "T", t, "Q", 1, f) for t in temps]) / 1e3
    ps = np.array([PropsSI("P", "T", t, "Q", 0, f) for t in temps]) / 1e5
    pg = np.array([PropsSI("P", "T", t, "Q", 1, f) for t in temps]) / 1e5  # differs from ps for blends
    pts = [result.states[k] for k in ("1", "2", "3", "4", "1")]
    h = [p.enthalpy_j_kg / 1e3 for p in pts]
    pressure = [p.pressure_pa / 1e5 for p in pts]
    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.plot(hf, ps, "k-", lw=1.2, label="saturated liquid")
    ax.plot(hg, pg, "k--", lw=1.2, label="saturated vapor")
    ax.plot(h, pressure, "o-", lw=2.2, color="#0072B2", label="cycle")
    for key in ("1", "2", "3", "4"):
        p = result.states[key]
        ax.annotate(key, (p.enthalpy_j_kg / 1e3, p.pressure_pa / 1e5), xytext=(6, 5), textcoords="offset points")
    ax.set_yscale("log")
    # Frame the cycle rather than the whole dome.
    ax.set_ylim(min(pressure) / 2.5, max(pressure) * 2.5)
    ax.set_xlim(min(h) - 60, max(h) + 40)
    ax.set_xlabel("Specific enthalpy (kJ/kg)")
    ax.set_ylabel("Pressure (bar, log scale)")
    ax.set_title(f"{f} vapor-compression cycle")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_hx_profiles(result: CoupledResult, path: str | Path) -> None:
    fig = hx_profiles_figure(result)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def hx_profiles_figure(result):
    """Return a Figure with evaporator and condenser temperature profiles.

    Accepts any result with ``evaporator`` and ``condenser`` moving-boundary results.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, hx, title in zip(axes, (result.evaporator, result.condenser), ("Evaporator", "Condenser")):
        x = np.linspace(0, 1, len(hx.refrigerant_temperature_c))
        ax.plot(x, hx.refrigerant_temperature_c, lw=2, label="refrigerant")
        ax.plot(x, hx.secondary_temperature_c, lw=2, label="secondary fluid")
        ax.set_title(f"{title}: min approach {hx.minimum_approach_k:.2f} K")
        ax.set_xlabel("Position along refrigerant flow")
        ax.set_ylabel("Temperature (°C)")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.tight_layout()
    return fig


def plot_cop_sweep(base: CycleInputs, path: str | Path) -> None:
    condensing = np.linspace(max(base.evaporating_temperature_c + 12, 25), 60, 36)
    cop = []
    carnot = []
    for tc in condensing:
        r = solve_cycle(CycleInputs(**{**base.__dict__, "condensing_temperature_c": float(tc)}))
        cop.append(r.cop_cooling)
        carnot.append(r.carnot_cop_cooling)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(condensing, cop, lw=2.2, label="modeled COP")
    ax.plot(condensing, carnot, "--", lw=1.5, label="Carnot upper bound")
    ax.set_xlabel("Condensing saturation temperature (°C)")
    ax.set_ylabel("Cooling COP")
    ax.set_title("Effect of temperature lift on COP")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)

