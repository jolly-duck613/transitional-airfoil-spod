"""Project-specific input/output helpers.

Keep solver-specific parsing in this module so that numerical analysis code does
not depend on filenames or a personal directory layout. Never embed credentials
or proprietary data in source files.
"""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray


def load_snapshots(path: str | Path) -> tuple[NDArray[np.float64], dict]:
    """Load time-resolved CFD states and their metadata.

    Parameters
    ----------
    path:
        Approved local data file or directory. Do not point public examples at
        private laboratory storage.

    Returns
    -------
    states:
        Array shaped ``(n_time, n_space)``. Flatten multiple variables in a
        documented, consistent order.
    metadata:
        Grid shape, coordinates, variable names, units, sampling interval,
        flattening convention, masks, and any spatial quadrature weights.

    Notes
    -----
    Replace this placeholder with a reader for your authorised data format.
    Validate uniform sampling and consistent meshes before spectral analysis.
    """
    data_path = Path(path)
    raise NotImplementedError(
        f"Add an approved data reader for {data_path}. "
        "See sample_data/README.md for the expected contract."
    )


def save_analysis(path: str | Path, **arrays: NDArray) -> None:
    """Save small derived arrays; avoid writing raw or restricted snapshots."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **arrays)
