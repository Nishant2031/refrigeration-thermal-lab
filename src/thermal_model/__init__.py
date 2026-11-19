"""Educational vapor-compression and heat-exchanger models (SI units)."""

from .cycle import CycleInputs, CycleResult, solve_cycle
from .coupled import CoupledInputs, CoupledResult, solve_coupled_system
from .heat_exchangers import (
    PhaseChangeResult,
    SensibleHXResult,
    effectiveness_ntu,
    lmtd,
    phase_change_hx,
    sensible_hx,
)
from .evaporator import (LiquidTaskResult, SaturationScreeningResult,
    RefrigerantFixtureResult, liquid_cooling_task,
    constant_saturation_evaporator, refrigerant_evaporator_states)
from .map_coupled import MapCoupledInputs, MapCoupledResult, solve_map_coupled_system
from .plate_evaporator import PlateEvaporatorGeometry, PlateEvaporatorResult, plate_evaporator_screening

__all__ = [
    "CycleInputs", "CycleResult", "solve_cycle",
    "CoupledInputs", "CoupledResult", "solve_coupled_system",
    "PhaseChangeResult", "SensibleHXResult", "effectiveness_ntu",
    "lmtd", "phase_change_hx", "sensible_hx",
    "LiquidTaskResult", "SaturationScreeningResult", "RefrigerantFixtureResult",
    "liquid_cooling_task", "constant_saturation_evaporator", "refrigerant_evaporator_states",
    "MapCoupledInputs", "MapCoupledResult", "solve_map_coupled_system",
    "PlateEvaporatorGeometry", "PlateEvaporatorResult", "plate_evaporator_screening",
]

