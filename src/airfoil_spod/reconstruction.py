"""Helpers for interpreting complex SPOD modes and coherent reconstructions."""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def mode_at_phase(mode: ArrayLike, phase_radians: float) -> NDArray[np.float64]:
    """Return ``Re(phi exp(i theta))`` for a complex spatial mode.

    A mode's absolute phase is arbitrary; use phase sequences to show relative
    evolution rather than assigning physical meaning to a single phase origin.
    """
    complex_mode = np.asarray(mode, dtype=complex)
    return np.real(complex_mode * np.exp(1j * phase_radians))


def phase_sequence(mode: ArrayLike, n_phases: int = 8) -> NDArray[np.float64]:
    """Create equally spaced views across one oscillation cycle."""
    if n_phases < 2:
        raise ValueError("n_phases must be at least two")
    phases = np.linspace(0.0, 2.0 * np.pi, n_phases, endpoint=False)
    return np.stack([mode_at_phase(mode, phase) for phase in phases])


def harmonic_reconstruction(
    mode: ArrayLike,
    *,
    coefficient: complex,
    frequency: float,
    times: ArrayLike,
) -> NDArray[np.float64]:
    """Reconstruct one idealised frequency-mode contribution over time.

    This is not a complete inverse-SPOD reconstruction. A research reconstruction
    must use coefficients obtained with a documented projection or oblique-
    projection method and must account for all selected blocks/frequencies.
    """
    complex_mode = np.asarray(mode, dtype=complex)
    time_array = np.asarray(times, dtype=float)
    oscillation = coefficient * np.exp(2j * np.pi * frequency * time_array)
    return np.real(oscillation[:, None] * complex_mode[None, :])
