"""Plotting helpers for PSD and SPOD spectra."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike


def plot_spectrum(
    frequencies: ArrayLike,
    eigenvalues: ArrayLike,
    output: str | Path,
    *,
    max_modes: int = 3,
) -> None:
    """Save leading SPOD eigenvalue curves without assuming a flow normalisation."""
    f = np.asarray(frequencies)
    values = np.asarray(eigenvalues)
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(7, 4))
    for index in range(min(max_modes, values.shape[1])):
        axis.semilogy(f, np.maximum(values[:, index], np.finfo(float).tiny), label=f"Mode {index + 1}")
    axis.set(xlabel="Frequency", ylabel="SPOD eigenvalue", title="SPOD spectrum")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_psd(frequencies: ArrayLike, psd: ArrayLike, output: str | Path) -> None:
    """Save a Welch PSD used as one frequency-domain validation check."""
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(7, 4))
    axis.semilogy(frequencies, np.maximum(psd, np.finfo(float).tiny))
    axis.set(xlabel="Frequency", ylabel="PSD", title="Welch power spectral density")
    axis.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
