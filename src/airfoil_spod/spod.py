"""Compact snapshot-form SPOD implementation for modest educational examples.

For a production CFD workflow, verify normalisation and spatial weighting against
a trusted implementation and use a memory-aware method suitable for the dataset.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import get_window

from .preprocessing import validate_snapshot_matrix


@dataclass(frozen=True)
class SPODResult:
    """Frequency axis, modal energies, complex modes, and analysis metadata.

    ``eigenvalues`` has shape ``(n_frequency, n_mode)`` and ``modes`` has shape
    ``(n_frequency, n_mode, n_space)``.
    """

    frequencies: NDArray[np.float64]
    eigenvalues: NDArray[np.float64]
    modes: NDArray[np.complex128]
    block_size: int
    overlap: int
    n_blocks: int


def _windowed_fft_blocks(
    states: NDArray[np.float64],
    *,
    block_size: int,
    overlap: int,
    window: str,
) -> tuple[NDArray[np.complex128], float]:
    """Return one-sided FFT blocks shaped ``(frequency, block, space)``."""
    step = block_size - overlap
    starts = np.arange(0, states.shape[0] - block_size + 1, step)
    if starts.size < 2:
        raise ValueError("SPOD requires at least two complete Welch blocks")

    taper = get_window(window, block_size, fftbins=True).astype(float)
    window_energy = float(np.sum(taper**2))
    blocks = []
    for start in starts:
        block = states[start : start + block_size]
        block = block - block.mean(axis=0, keepdims=True)
        blocks.append(np.fft.rfft(block * taper[:, None], axis=0))
    return np.stack(blocks, axis=1), window_energy


def compute_spod(
    states: ArrayLike,
    *,
    dt: float,
    block_size: int,
    overlap: int,
    n_modes: int | None = None,
    weights: ArrayLike | None = None,
    window: str = "hann",
) -> SPODResult:
    """Compute SPOD modes using the method of snapshots.

    Parameters
    ----------
    states:
        Mean-subtracted or raw states shaped ``(n_time, n_space)``. Each block is
        detrended by removing its temporal mean.
    dt:
        Uniform time between snapshots.
    block_size, overlap, window:
        Welch segmentation settings.
    n_modes:
        Number of leading modes retained at each frequency. It cannot exceed the
        number of blocks.
    weights:
        Positive spatial quadrature/energy weights. Unit weights are only
        appropriate for equally weighted degrees of freedom.

    Notes
    -----
    Eigenvalues use a consistent density-like scaling for this scaffold. Establish
    and test the exact convention required for comparison with another library.
    """
    array = validate_snapshot_matrix(states)
    if dt <= 0:
        raise ValueError("dt must be positive")
    if not 0 <= overlap < block_size <= array.shape[0]:
        raise ValueError("require 0 <= overlap < block_size <= n_time")

    fft_blocks, window_energy = _windowed_fft_blocks(
        array, block_size=block_size, overlap=overlap, window=window
    )
    n_frequencies, n_blocks, n_space = fft_blocks.shape
    retained = n_blocks if n_modes is None else min(n_modes, n_blocks)
    if retained < 1:
        raise ValueError("n_modes must be positive")

    spatial_weights = (
        np.ones(n_space, dtype=float)
        if weights is None
        else np.asarray(weights, dtype=float)
    )
    if spatial_weights.shape != (n_space,) or np.any(spatial_weights <= 0):
        raise ValueError("weights must contain one positive value per spatial degree of freedom")

    eigenvalues = np.zeros((n_frequencies, retained), dtype=float)
    modes = np.zeros((n_frequencies, retained, n_space), dtype=complex)
    scale = dt / (window_energy * n_blocks)

    for frequency_index in range(n_frequencies):
        q_hat = fft_blocks[frequency_index].T  # (space, blocks)
        correlation = scale * (q_hat.conj().T @ (spatial_weights[:, None] * q_hat))
        values, vectors = np.linalg.eigh(correlation)
        order = np.argsort(values)[::-1]
        values = np.maximum(values[order][:retained].real, 0.0)
        vectors = vectors[:, order][:, :retained]

        positive = values > np.finfo(float).eps
        frequency_modes = np.zeros((n_space, retained), dtype=complex)
        frequency_modes[:, positive] = (
            np.sqrt(scale)
            * q_hat
            @ vectors[:, positive]
            / np.sqrt(values[positive])[None, :]
        )
        eigenvalues[frequency_index] = values
        modes[frequency_index] = frequency_modes.T

    frequencies = np.fft.rfftfreq(block_size, d=dt)
    return SPODResult(
        frequencies=frequencies,
        eigenvalues=eigenvalues,
        modes=modes,
        block_size=block_size,
        overlap=overlap,
        n_blocks=n_blocks,
    )
