"""CoolProp-backed thermodynamic state helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite

from CoolProp.CoolProp import PhaseSI, PropsSI


@dataclass(frozen=True)
class StatePoint:
    """A refrigerant state. All stored quantities use base SI units."""

    label: str
    pressure_pa: float
    temperature_k: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float
    density_kg_m3: float
    quality: float | None
    phase: str

    @classmethod
    def from_pair(
        cls, label: str, fluid: str, key1: str, value1: float,
        key2: str, value2: float,
    ) -> "StatePoint":
        p = float(PropsSI("P", key1, value1, key2, value2, fluid))
        t = float(PropsSI("T", key1, value1, key2, value2, fluid))
        h = float(PropsSI("Hmass", key1, value1, key2, value2, fluid))
        s = float(PropsSI("Smass", key1, value1, key2, value2, fluid))
        rho = float(PropsSI("Dmass", key1, value1, key2, value2, fluid))
        phase = str(PhaseSI(key1, value1, key2, value2, fluid))
        try:
            raw_q = float(PropsSI("Q", key1, value1, key2, value2, fluid))
            q = raw_q if isfinite(raw_q) and 0.0 <= raw_q <= 1.0 else None
        except ValueError:
            q = None
        return cls(label, p, t, h, s, rho, q, phase)

    def as_display_dict(self) -> dict[str, float | str | None]:
        return {
            "label": self.label,
            "pressure_kPa": self.pressure_pa / 1e3,
            "temperature_C": self.temperature_k - 273.15,
            "enthalpy_kJ_kg": self.enthalpy_j_kg / 1e3,
            "entropy_kJ_kgK": self.entropy_j_kg_k / 1e3,
            "density_kg_m3": self.density_kg_m3,
            "quality": self.quality,
            "phase": self.phase,
        }

    def as_si_dict(self) -> dict:
        return asdict(self)

