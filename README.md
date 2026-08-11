# SPOD Analysis of Transitional Airfoil Flow

A compact Python repository for Fast Fourier Transform (FFT), Welch power spectral density (PSD), and Spectral Proper Orthogonal Decomposition (SPOD) analysis.

The repository is self-contained: it includes the core numerical routines, plotting helpers, a synthetic smoke test, and automated tests. It does not contain CFD datasets or research results.

## Included functionality

- snapshot validation and temporal mean subtraction;
- Welch PSD estimation;
- chord-based Strouhal-number conversion;
- snapshot-form SPOD with optional spatial weights;
- phase views and simple harmonic reconstruction;
- PSD and SPOD-spectrum plotting;
- a reproducible synthetic-data smoke test.

## Repository structure

```text
transitional-airfoil-spod/
├── README.md
├── LICENSE
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── src/
│   └── airfoil_spod/
│       ├── __init__.py
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

## Verify the installation

Run the automated tests:

```bash
python -m unittest discover -s tests -v
```

Run the self-contained smoke test:

```bash
python scripts/run_example.py --output results/example
```

It creates synthetic travelling waves, computes their Welch PSD and SPOD spectrum, and writes two plots to `results/example/`.

## Basic use

The numerical routines expect a real-valued snapshot matrix shaped `(n_time, n_space)`.

```python
from airfoil_spod.preprocessing import subtract_mean
from airfoil_spod.spectral import welch_psd
from airfoil_spod.spod import compute_spod

fluctuations, mean = subtract_mean(snapshot_matrix)

frequencies, psd = welch_psd(
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

For nonuniform meshes, pass positive spatial quadrature weights to `compute_spod`. Choose the sampling interval, block size, overlap, and window consistently with the analysis.

## Frequency convention

Chord-based Strouhal number is calculated as

```text
St_c = f c / U_inf
```

where `f` is frequency, `c` is chord, and `U_inf` is freestream velocity.

## Licence

Released under the MIT License. Copyright © 2026 Dov Noimark.

