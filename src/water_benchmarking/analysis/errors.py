"""Error estimates for averages over correlated time series.

Consecutive MD frames are not independent, so the naive standard error is far too
small.  Block averaging inflates the block size until the block-mean variance
stops growing, which is where the blocks have become effectively independent.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Estimate:
    """A mean with an honest uncertainty."""

    mean: float
    error: float
    n_samples: int

    def __str__(self) -> str:
        return f"{self.mean:.4g} +/- {self.error:.2g}"


def block_average(values: np.ndarray, min_blocks: int = 8) -> Estimate:
    """Mean and standard error of a correlated series by block averaging.

    The error is taken from the plateau of the block-size curve -- in practice the
    largest block count that still leaves min_blocks blocks, which is where the
    estimate has stopped rising but has not yet gone noisy.
    """
    values = np.asarray(values, dtype=float).ravel()
    n = values.size
    if n < min_blocks:
        return Estimate(float(values.mean()), float("nan"), n)

    errors = []
    size = 1
    while n // size >= min_blocks:
        n_blocks = n // size
        blocks = values[: n_blocks * size].reshape(n_blocks, size).mean(axis=1)
        errors.append(blocks.std(ddof=1) / np.sqrt(n_blocks))
        size *= 2
    return Estimate(float(values.mean()), float(max(errors)), n)


def drop_equilibration(values: np.ndarray, fraction: float = 0.1) -> np.ndarray:
    """Discard the leading fraction of a series.

    The run ladder already equilibrates before production, so this only guards
    against a residual drift in the first segment.
    """
    values = np.asarray(values, dtype=float)
    return values[int(len(values) * fraction):]
