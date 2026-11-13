"""Owned immutable preparation and explicit Gaussian survey loop closure.

Coordinates are local ENU metres. Buffers are copied before solving; callers
must not mutate exporters concurrently with preparation.
"""

from __future__ import annotations

import struct
import sys
from enum import StrEnum
from types import MappingProxyType

from cave_closure._native import (
    CancellationToken,
    IncompleteResultError,
    LoopClosureAlgorithm,
    Prepared,
    Report,
    SolveError,
    ValidationError,
    prepare_arrays,
    prepare_json,
)

__all__ = [
    "CancellationToken",
    "LoopClosureMode",
    "LoopClosureAlgorithm",
    "Prepared",
    "Report",
    "ValidationError",
    "SolveError",
    "IncompleteResultError",
    "coordinates",
    "coordinate_buffer",
    "prepare_arrays",
    "prepare_buffers",
    "prepare_json",
]


class LoopClosureMode(StrEnum):
    """Loop-closure modes accepted by the array and buffer preparation APIs."""

    SPATIAL = "spatial"
    HORIZONTAL = "horizontal"


def coordinates(result: Report):
    """Immutable coordinate mapping; values are immutable triples."""
    return MappingProxyType({key: tuple(value) for key, value in result.coordinates})


def _copy_buffer(value, *, name, shape, code, itemsize):
    view = memoryview(value)
    if not view.c_contiguous:
        raise ValueError(f"{name}: expected a C-contiguous buffer")
    if view.shape != shape or view.itemsize != itemsize:
        raise ValueError(f"{name}: invalid shape or item size")
    # Accept only unprefixed/native or explicitly matching native endian formats.
    native_endian = "<" if sys.byteorder == "little" else ">"
    codes = {code, "L"} if code == "Q" else {code}
    allowed = {
        prefix + item for item in codes for prefix in ("", "@", "=", native_endian)
    }
    if view.format not in allowed:
        raise ValueError(f"{name}: unsupported format {view.format!r}")
    # bytes() performs the owned copy; memoryview.cast avoids exporter-specific
    # indexing and handles accepted explicit endian spellings consistently.
    copied = bytes(view)
    flat = memoryview(copied).cast(code).tolist()
    if len(shape) == 2:
        return [flat[i : i + shape[1]] for i in range(0, len(flat), shape[1])]
    width = shape[1] * shape[2]
    return [
        [flat[i + j : i + j + shape[2]] for j in range(0, width, shape[2])]
        for i in range(0, len(flat), width)
    ]


def prepare_buffers(
    *,
    station_ids,
    observation_ids,
    endpoints,
    displacements,
    covariances,
    mode: LoopClosureMode,
    algorithm: LoopClosureAlgorithm = LoopClosureAlgorithm.LSE_REDUCED_PCG,
    preserved_z=None,
    controls=(),
    relationships=(),
    reading_ids=None,
    provenance=(),
    require_control=False,
    allow_shared_reading_approximation=False,
):
    """Copy C-contiguous native uint64 endpoints and float64 vectors/covariances.

    Shapes: (observations, 2), (observations, 3), (observations, 3, 3).
    The full spatial covariance is supplied even in marginal-horizontal mode.
    NumPy is optional. Mutable exporters must remain unchanged during this call.
    """
    station_ids = tuple(station_ids)
    observation_ids = tuple(observation_ids)
    if len(station_ids) > 1_000_000 or len(observation_ids) > 4_000_000:
        raise ValueError("configured input resource limit exceeded")
    m = len(observation_ids)
    endpoints = _copy_buffer(
        endpoints, name="endpoints", shape=(m, 2), code="Q", itemsize=8
    )
    displacements = _copy_buffer(
        displacements, name="displacements", shape=(m, 3), code="d", itemsize=8
    )
    covariances = _copy_buffer(
        covariances, name="covariances", shape=(m, 3, 3), code="d", itemsize=8
    )
    return prepare_arrays(
        station_ids,
        observation_ids,
        endpoints,
        displacements,
        covariances,
        mode=mode,
        algorithm=algorithm,
        preserved_z=preserved_z,
        controls=list(controls),
        relationships=list(relationships),
        reading_ids=reading_ids,
        provenance=list(provenance),
        require_control=require_control,
        allow_shared_reading_approximation=allow_shared_reading_approximation,
    )


def coordinate_buffer(result: Report):
    """Read-only owned copy: three native float64 values per sorted station.

    The memoryview retains its bytes owner independently of the result lifetime.
    This convenience export makes no zero-copy guarantee.
    """
    values = [value for _, point in result.coordinates for value in point]
    return memoryview(struct.pack(f"{len(values)}d", *values)).cast("d")
