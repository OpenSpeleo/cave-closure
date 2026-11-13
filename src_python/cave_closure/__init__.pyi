from collections.abc import Mapping, Sequence
from enum import StrEnum
from typing import Any

from ._native import CancellationToken as CancellationToken
from ._native import IncompleteResultError as IncompleteResultError
from ._native import LoopClosureAlgorithm as LoopClosureAlgorithm
from ._native import Prepared as Prepared
from ._native import Report as Report
from ._native import SolveError as SolveError
from ._native import ValidationError as ValidationError
from ._native import prepare_arrays as prepare_arrays
from ._native import prepare_json as prepare_json

class Mode(StrEnum):
    SPATIAL = "spatial"
    HORIZONTAL = "horizontal"

def coordinates(result: Report) -> Mapping[str, tuple[float, float, float]]: ...
def prepare_buffers(
    *,
    station_ids: Sequence[str],
    observation_ids: Sequence[str],
    endpoints: Any,
    displacements: Any,
    covariances: Any,
    mode: Mode,
    algorithm: LoopClosureAlgorithm = ...,
    preserved_z: Sequence[float] | None = ...,
    controls: Sequence[Any] = ...,
    relationships: Sequence[Any] = ...,
    reading_ids: Sequence[Sequence[str]] | None = ...,
    provenance: Sequence[str] = ...,
    require_control: bool = ...,
    allow_shared_reading_approximation: bool = ...,
) -> Prepared: ...
def coordinate_buffer(result: Report) -> memoryview: ...
