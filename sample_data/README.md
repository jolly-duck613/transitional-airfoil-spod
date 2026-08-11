# Sample data

No research or laboratory CFD data is distributed with this repository.

The example script generates synthetic data in memory. Replace the placeholder loader in `src/airfoil_spod/io.py` only after confirming that the chosen input data may be used and shared.

## Suggested public example format

If you have permission to publish a small example, prefer a compact, clearly documented `.npz` file containing only the fields necessary to demonstrate the workflow:

```text
sample_snapshot_series.npz
├── states      # shape: (n_time, n_space), or documented equivalent
├── time        # shape: (n_time,)
├── x           # spatial coordinates
├── y           # spatial coordinates
└── metadata    # variable names, units, flattening order, and provenance
```

Document:

- whether values are dimensional or nondimensional;
- the sampling interval and whether sampling is uniform;
- the Reynolds and Mach numbers, angle of attack, and reference scales, if approved;
- which variables are included (`u`, `v`, density, vorticity, and so on);
- grid topology, masks, quadrature weights, and flattening convention;
- preprocessing already applied;
- the dataset licence and attribution requirements.

## Do not commit

- full solver outputs or raw simulation campaigns;
- unpublished geometries or operating conditions;
- data owned by a laboratory, sponsor, collaborator, or institution without approval;
- personal paths, credentials, access tokens, or machine configuration;
- files whose licence or provenance is uncertain.

For large approved datasets, publish them in an appropriate data repository and provide a DOI or download script instead of storing them in Git.
