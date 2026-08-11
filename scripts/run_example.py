"""Run a small, non-proprietary SPOD demonstration on synthetic travelling waves."""

import argparse
from pathlib import Path

import numpy as np

from airfoil_spod.plotting import plot_psd, plot_spectrum
from airfoil_spod.preprocessing import subtract_mean
from airfoil_spod.spectral import welch_psd
from airfoil_spod.spod import compute_spod


def synthetic_states(
    *, n_time: int = 1024, n_space: int = 80, dt: float = 0.02, seed: int = 7
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create two travelling waves plus noise for a reproducible smoke test."""
    rng = np.random.default_rng(seed)
    time = np.arange(n_time) * dt
    x = np.linspace(0.0, 1.0, n_space, endpoint=False)
    primary = np.sin(2 * np.pi * (3.0 * time[:, None] - 2.0 * x[None, :]))
    secondary = 0.35 * np.sin(2 * np.pi * (7.0 * time[:, None] - 5.0 * x[None, :]) + 0.4)
    noise = 0.12 * rng.standard_normal((n_time, n_space))
    return primary + secondary + noise, time, x


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/example"))
    args = parser.parse_args()

    dt = 0.02
    states, _, _ = synthetic_states(dt=dt)
    fluctuations, _ = subtract_mean(states)
    frequencies, density = welch_psd(
        fluctuations[:, 10], dt=dt, block_size=256, overlap=128
    )
    result = compute_spod(
        fluctuations, dt=dt, block_size=256, overlap=128, n_modes=3
    )

    args.output.mkdir(parents=True, exist_ok=True)
    plot_psd(frequencies, density, args.output / "welch_psd.png")
    plot_spectrum(result.frequencies, result.eigenvalues, args.output / "spod_spectrum.png")

    dominant_index = 1 + int(np.argmax(result.eigenvalues[1:, 0]))
    print(f"Wrote example figures to {args.output}")
    print(f"Leading non-zero-frequency SPOD peak: {result.frequencies[dominant_index]:.3f}")


if __name__ == "__main__":
    main()
