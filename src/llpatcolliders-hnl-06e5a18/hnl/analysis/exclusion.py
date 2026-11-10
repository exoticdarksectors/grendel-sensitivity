"""
analysis/exclusion.py

Extract (m_N, U²) exclusion regions from N_signal scans.

When both threshold crossings are resolved, the exclusion region is an
"island" in (m_N, U²) space:
- Too small U²: insufficient production rate
- Too large U²: HNL decays before reaching the detector

Either edge can remain open when the scan range does not reach the crossing.
Open edges are represented by NaN boundaries plus explicit boolean flags; the
scan endpoint must not be presented as a physical exclusion boundary.
"""

import numpy as np

from analysis.constants import N_THRESHOLD


def find_exclusion_band(u2_grid, N_grid, N_threshold=N_THRESHOLD):
    """
    Find the U² range where N_signal >= N_threshold.

    Parameters
    ----------
    u2_grid : (M,) log-spaced U² values
    N_grid : (M,) corresponding N_signal values
    N_threshold : float — minimum signal events for exclusion

    Returns
    -------
    dict with:
        u2_min : float — resolved lower edge (NaN if absent or open)
        u2_max : float — resolved upper edge (NaN if absent or open)
        u2_min_open : bool — exclusion continues below the scan range
        u2_max_open : bool — exclusion continues above the scan range
        peak_N : float — maximum N_signal over the scan
        peak_u2 : float — U² at peak sensitivity
        has_sensitivity : bool
    """
    peak_idx = np.argmax(N_grid)
    peak_N = float(N_grid[peak_idx])
    peak_u2 = float(u2_grid[peak_idx])

    mask = N_grid >= N_threshold
    if not np.any(mask):
        return {
            "u2_min": np.nan, "u2_max": np.nan,
            "u2_min_open": False, "u2_max_open": False,
            "peak_N": peak_N, "peak_u2": peak_u2,
            "has_sensitivity": False,
        }

    idx = np.where(mask)[0]
    i_lo, i_hi = idx[0], idx[-1]

    u2_min_open = i_lo == 0
    u2_max_open = i_hi == len(u2_grid) - 1
    u2_min = (
        np.nan if u2_min_open
        else _interpolate_threshold(u2_grid, N_grid, i_lo, N_threshold, "lower")
    )
    u2_max = (
        np.nan if u2_max_open
        else _interpolate_threshold(u2_grid, N_grid, i_hi, N_threshold, "upper")
    )

    return {
        "u2_min": u2_min, "u2_max": u2_max,
        "u2_min_open": u2_min_open, "u2_max_open": u2_max_open,
        "peak_N": peak_N, "peak_u2": peak_u2,
        "has_sensitivity": True,
    }


def _interpolate_threshold(u2, N, idx, N_thr, side):
    """Interpolate in log(U²) space to find where N crosses the threshold."""
    if side == "lower" and idx > 0:
        N_below, N_above = N[idx - 1], N[idx]
        if N_above > N_below:
            frac = (N_thr - N_below) / (N_above - N_below)
            log_lo = np.log10(u2[idx - 1])
            log_hi = np.log10(u2[idx])
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    elif side == "upper" and idx < len(u2) - 1:
        N_above, N_below = N[idx], N[idx + 1]
        if N_above > N_below:
            frac = (N_above - N_thr) / (N_above - N_below)
            log_lo = np.log10(u2[idx])
            log_hi = np.log10(u2[idx + 1])
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    return float(u2[idx])
