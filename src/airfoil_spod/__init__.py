"""Tools for a compact, educational SPOD analysis workflow."""

from .preprocessing import subtract_mean
from .spectral import welch_psd
from .spod import SPODResult, compute_spod

__all__ = ["SPODResult", "compute_spod", "subtract_mean", "welch_psd"]
