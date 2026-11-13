from collections.abc import Mapping, Sequence
from typing import Any

class ValidationError(ValueError): ...
class SolveError(RuntimeError): ...
class IncompleteResultError(SolveError): ...

class CancellationToken:
    def __init__(self) -> None: ...
    def cancel(self) -> None: ...
    @property
    def cancelled(self) -> bool: ...

class Prepared:
    def solve(
        self,
        *,
        max_iterations: int = ...,
        reduce: bool = ...,
        cancellation: CancellationToken | None = ...,
    ) -> Report: ...

class Report:
    @property
    def status(self) -> str: ...
    @property
    def coordinates(self) -> list[tuple[str, list[float]]]: ...
    @property
    def observations(self) -> list[tuple[str, list[float]]]: ...
    def require_converged(self) -> None: ...
    def to_json(self) -> str: ...

def coordinates(result: Report) -> Mapping[str, tuple[float, float, float]]: ...
def prepare_json(
    input: str,
    *,
    require_control: bool = ...,
    allow_shared_reading_approximation: bool = ...,
    max_stations: int = ...,
    max_observations: int = ...,
) -> Prepared: ...
def prepare_arrays(
    station_ids: Sequence[str],
    observation_ids: Sequence[str],
    endpoints: Sequence[Sequence[int]],
    displacements: Sequence[Sequence[float]],
    covariances: Sequence[Sequence[Sequence[float]]],
    *,
    mode: str,
    preserved_z: Sequence[float] | None = ...,
    controls: Sequence[Any] = ...,
    relationships: Sequence[Any] = ...,
    reading_ids: Sequence[Sequence[str]] | None = ...,
    provenance: Sequence[str] = ...,
    require_control: bool = ...,
    allow_shared_reading_approximation: bool = ...,
) -> Prepared: ...
def prepare_buffers(
    *,
    station_ids: Sequence[str],
    observation_ids: Sequence[str],
    endpoints: Any,
    displacements: Any,
    covariances: Any,
    mode: str,
    preserved_z: Sequence[float] | None = ...,
    controls: Sequence[Any] = ...,
    relationships: Sequence[Any] = ...,
    reading_ids: Sequence[Sequence[str]] | None = ...,
    provenance: Sequence[str] = ...,
    require_control: bool = ...,
    allow_shared_reading_approximation: bool = ...,
) -> Prepared: ...
def coordinate_buffer(result: Report) -> memoryview: ...
