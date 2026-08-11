"""Welch spectral estimates used to validate SPOD peak frequencies."""

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import welch


def welch_psd(
    signal: ArrayLike,
    *,
    dt: float,
    block_size: int,
    overlap: int,
    window: str = "hann",
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Estimate a one-sided PSD using Welch averaging.

    The returned density uses cycles per unit time. Compare PSD peaks with SPOD
    peaks using identical sampling and compatible block/window choices, while
    remembering that the two estimates answer related but different questions.
    """
    values = np.asarray(signal, dtype=float)
    if values.ndim != 1:
        raise ValueError("signal must be one-dimensional")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if not 0 <= overlap < block_size <= values.size:
        raise ValueError("require 0 <= overlap < block_size <= signal length")

    frequencies, density = welch(
        values,
        fs=1.0 / dt,
        window=window,
        nperseg=block_size,
        noverlap=overlap,
        detrend="constant",
        scaling="density",
        return_onesided=True,
    )
    return frequencies, density


def strouhal_number(
    frequencies: ArrayLike, *, chord: float, freestream_velocity: float
) -> NDArray[np.float64]:
    """Convert frequencies to chord-based Strouhal numbers ``St = f c / U_inf``."""
    if chord <= 0 or freestream_velocity <= 0:
        raise ValueError("chord and freestream_velocity must be positive")
    return np.asarray(frequencies, dtype=float) * chord / freestream_velocity
