"""Functional redundancy and core-function classification.

Data required
-------------
abundance : samples x MAG abundance matrix
carriers  : Boolean vector identifying MAGs carrying one function
"""

import numpy as np


def effective_mag_number(x):
    """Inverse Simpson effective number of contributing MAGs."""
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if x.size == 0:
        return np.nan
    p = x / x.sum()
    return 1.0 / np.sum(p ** 2)


def summarize_function(abundance, carriers, core_threshold=0.01):
    x = np.asarray(abundance, dtype=float)[:, np.asarray(carriers, dtype=bool)]
    summed = x.sum(axis=1)
    effective = np.array([effective_mag_number(row) for row in x])

    return {
        "min_abundance": float(np.min(summed)),
        "median_abundance": float(np.median(summed)),
        "max_abundance": float(np.max(summed)),
        "samples_at_or_above_threshold": int(np.sum(summed >= core_threshold)),
        "core_all_samples": bool(np.all(summed >= core_threshold)),
        "median_effective_mag_number": float(np.nanmedian(effective)),
    }
