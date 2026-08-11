# SPOD Analysis of Transitional Airfoil Flow

An editable Python scaffold for analysing coherent structures in time-resolved computational fluid dynamics (CFD) data using Fast Fourier Transforms (FFT), Welch power spectral density (PSD) estimation, and Spectral Proper Orthogonal Decomposition (SPOD).

The project is designed around transitional, low-Reynolds-number airfoil flow, including laminar separation, recirculation, shear-layer transition, vortex shedding, and wake dynamics. It contains no proprietary CFD data, unpublished results, or solver-specific files.

> **Publication note:** Before publishing research code, figures, geometry, flow conditions, or results, confirm with your supervisor and institution that they may be shared.

## What this repository does

The intended analysis pipeline is:

1. Load a uniformly sampled sequence of CFD snapshots.
2. Assemble one or more flow variables into a snapshot matrix.
3. Remove the temporal mean to obtain fluctuations.
4. Divide the time series into overlapping, windowed Welch blocks.
5. Apply an FFT to every block.
6. Estimate conventional PSDs for frequency validation.
7. Solve the frequency-by-frequency SPOD eigenvalue problem.
8. Rank modes by spectral energy and identify dominant frequencies.
9. Visualise the leading modes at several phases.
10. Reconstruct selected coherent structures.
11. Compare frequencies, modes, and reconstructions with the original instantaneous flow fields.

## Repository structure

```text
transitional-airfoil-spod/
├── README.md
├── LICENSE
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── sample_data/
│   └── README.md
├── paper/
│   └── README.md
├── figures/
│   └── README.md
├── src/
│   └── airfoil_spod/
│       ├── __init__.py
│       ├── io.py
│       ├── preprocessing.py
│       ├── spectral.py
│       ├── spod.py
│       ├── reconstruction.py
│       └── plotting.py
├── scripts/
│   └── run_example.py
└── tests/
    └── test_pipeline.py
```

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Run the automated numerical checks with `python -m unittest discover -s tests -v`.

## Quick start with synthetic data

Run the included example without downloading any research data:

```bash
python scripts/run_example.py --output results/example
```

The script creates a small synthetic travelling-wave dataset, computes a reference Welch PSD and SPOD spectrum, and writes plots to `results/example/`. It is a structural demonstration, not a physical model of an airfoil.

## Adapting the scaffold to CFD snapshots

Implement the project-specific loader in `src/airfoil_spod/io.py`. The analysis routines expect a real-valued array with shape:

```text
(n_time, n_space)
```

For multiple variables, flatten and concatenate the spatial fields consistently, for example:

```python
state = np.concatenate([u.ravel(), v.ravel()])
```

Keep a record of the grid shape, mask, variables, units, and flattening order so that modes can be mapped back to physical space. For nonuniform meshes or multiple variables, provide physically meaningful spatial quadrature weights to `compute_spod`.

Example analysis:

```python
from airfoil_spod.preprocessing import subtract_mean
from airfoil_spod.spectral import welch_psd
from airfoil_spod.spod import compute_spod

fluctuations, mean = subtract_mean(snapshot_matrix)

f_psd, psd = welch_psd(
    probe_signal,
    dt=0.001,
    block_size=256,
    overlap=128,
)

result = compute_spod(
    fluctuations,
    dt=0.001,
    block_size=256,
    overlap=128,
    n_modes=10,
)
```

Choose `block_size`, overlap, window, sampling interval, and the number of statistically independent blocks deliberately. These choices control frequency resolution, estimator variance, and convergence.

## Nondimensional frequency

If `f` is frequency, `c` is airfoil chord, and `U_inf` is freestream velocity, the chord-based Strouhal number is

```text
St = f c / U_inf
```

State the reference length and velocity whenever reporting a Strouhal number.

## Validation methodology

SPOD results should not be accepted solely because a leading eigenvalue is large. Suggested checks include:

- Compare dominant SPOD frequencies with Welch PSD peaks from physically relevant probes or integrated signals.
- Repeat the analysis with reasonable changes to block length, overlap, window, and number of snapshots.
- Check that leading eigenvalue separation and spatial structures are sufficiently converged.
- Inspect several phases of each complex leading mode, for example `real(phi * exp(i theta))`.
- Reconstruct selected frequency-mode contributions and compare them with instantaneous velocity, vorticity, density, turbulent kinetic energy, and Reynolds-stress fields where appropriate.
- Project or correlate snapshots with the candidate coherent structure to select representative physical events.
- Confirm that spatial structures align with plausible mechanisms in the separation bubble, shear layer, recirculation region, and wake.
- Check units, FFT normalisation, one-sided spectral scaling, frequency ordering, grid weighting, and energy conservation on a controlled test case.
- Distinguish robust coherent structures from broadband turbulence and numerical artefacts.

Agreement between two methods that share the same preprocessing is useful but not fully independent validation. Where possible, combine spectral, spatial, numerical-convergence, and physical checks.

## Expected research outputs

Typical public-safe outputs may include:

- an SPOD eigenspectrum or modal-energy spectrum;
- Welch PSD comparisons at selected probes;
- leading-mode magnitude and phase views;
- phase-resolved views of a coherent structure;
- low-order modal reconstructions;
- comparisons with authorised instantaneous or mean fields.

Only include outputs that have been approved for public release.

## Reproducibility checklist

- Record the sampling interval and snapshot count.
- Record the analysed variables and their units.
- Record mean-removal and detrending choices.
- Record the window, block length, overlap, FFT convention, and one-sided scaling.
- Record mesh weighting and variable normalisation.
- Save configuration values alongside every generated figure.
- Use a fixed random seed for synthetic tests.
- Document software versions and any solver-specific preprocessing.

## Limitations

This is an educational scaffold, not a validated SPOD library. The included implementation is intentionally compact and suited to modest examples. Large CFD datasets typically require streaming I/O, parallel FFTs, memory-aware algorithms, and careful treatment of spatial inner products.

## Licence and citation

The scaffold is released under the MIT License. Add your name and the correct year to `LICENSE` before publication. If this code supports a paper, add the final citation and DOI to `paper/README.md` once public.

## Acknowledgements

Add approved acknowledgements for supervisors, collaborators, the CFD laboratory, funding bodies, and computing facilities here.
