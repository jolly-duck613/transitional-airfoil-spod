"""Preprocessing utilities for snapshot matrices."""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def validate_snapshot_matrix(states: ArrayLike) -> NDArray[np.float64]:
    """Return a finite 2-D floating-point array shaped ``(time, space)``."""
    array = np.asarray(states, dtype=float)
    if array.ndim != 2:
        raise ValueError("states must have shape (n_time, n_space)")
    if array.shape[0] < 2 or array.shape[1] < 1:
        raise ValueError("states must contain at least two times and one degree of freedom")
    if not np.all(np.isfinite(array)):
        raise ValueError("states contain NaN or infinite values")
    return array


def subtract_mean(
    states: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Subtract the temporal mean from every spatial degree of freedom."""
    array = validate_snapshot_matrix(states)
    mean = array.mean(axis=0)
    return array - mean, mean


def standardise_variables(
    states: ArrayLike, scales: ArrayLike
) -> NDArray[np.float64]:
    """Scale variables before concatenation when their units or magnitudes differ.

    Choose scales from the physical energy inner product or documented reference
    quantities, not merely to make a plot look balanced.
    """
    array = validate_snapshot_matrix(states)
    scale_array = np.asarray(scales, dtype=float)
    if scale_array.shape != (array.shape[1],):
        raise ValueError("scales must contain one value per spatial degree of freedom")
    if np.any(scale_array <= 0):
        raise ValueError("all scales must be positive")
    return array / scale_array
