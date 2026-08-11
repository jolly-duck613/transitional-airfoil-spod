"""Numerical smoke tests for the public synthetic-data workflow."""

import unittest

import numpy as np

from airfoil_spod.preprocessing import subtract_mean
from airfoil_spod.reconstruction import phase_sequence
from airfoil_spod.spectral import strouhal_number, welch_psd
from airfoil_spod.spod import compute_spod


class PipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(42)
        self.dt = 0.01
        time = np.arange(1024) * self.dt
        x = np.linspace(0.0, 1.0, 32, endpoint=False)
        self.states = np.sin(
            2 * np.pi * (5.0 * time[:, None] - 2.0 * x[None, :])
        ) + 0.03 * rng.standard_normal((1024, 32))

    def test_synthetic_frequency_and_shapes(self) -> None:
        fluctuations, _ = subtract_mean(self.states)
        frequencies, density = welch_psd(
            fluctuations[:, 0], dt=self.dt, block_size=256, overlap=128
        )
        result = compute_spod(
            fluctuations,
            dt=self.dt,
            block_size=256,
            overlap=128,
            n_modes=3,
        )

        welch_peak = frequencies[1 + np.argmax(density[1:])]
        spod_peak = result.frequencies[1 + np.argmax(result.eigenvalues[1:, 0])]
        self.assertAlmostEqual(welch_peak, 5.0, delta=0.4)
        self.assertAlmostEqual(spod_peak, 5.0, delta=0.4)
        self.assertTrue(np.allclose(fluctuations.mean(axis=0), 0.0, atol=1e-12))
        self.assertEqual(phase_sequence(result.modes[13, 0], 4).shape, (4, 32))

    def test_strouhal_conversion(self) -> None:
        result = strouhal_number([5.0], chord=1.0, freestream_velocity=10.0)
        np.testing.assert_allclose(result, [0.5])


if __name__ == "__main__":
    unittest.main()
