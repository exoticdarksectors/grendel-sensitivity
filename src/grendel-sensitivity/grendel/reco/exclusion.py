"""Extract exclusion regions in the (mass, coupling) plane from N_signal scans.

When both threshold crossings are resolved, the exclusion region is an
"island" in (m_N, U²) space:
- Too small U²: insufficient production rate
- Too large U²: the particle decays before reaching the detector

Either edge can remain open when the scan range does not reach the crossing.
Open edges are represented by NaN boundaries plus explicit boolean flags; the
scan endpoint must not be presented as a physical exclusion boundary.
"""

import numpy as np

from ..constants import N_THRESHOLD


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


def find_exclusion_band_refined(
    u2_grid,
    N_grid,
    evaluate,
    N_threshold=N_THRESHOLD,
    log10_tolerance=1.0e-6,
):
    """Refine a coupling island using a deterministic yield evaluator.

    The coarse scan locates the peak and crossing brackets.  The peak is then
    maximized and each closed edge is bisected by evaluating the already-built
    acceptance Monte Carlo, without sampling any new events.
    """
    u2 = np.asarray(u2_grid, dtype=float)
    signal = np.asarray(N_grid, dtype=float)
    if u2.ndim != 1 or signal.shape != u2.shape or len(u2) < 3:
        raise ValueError("u2_grid and N_grid must be equal-length 1D arrays")
    if np.any(~np.isfinite(u2)) or np.any(u2 <= 0.0):
        raise ValueError("u2_grid must contain finite positive couplings")
    if np.any(np.diff(u2) <= 0.0):
        raise ValueError("u2_grid must be strictly increasing")

    cache = {float(x): float(y) for x, y in zip(u2, signal)}

    def value(x):
        x = float(x)
        if x not in cache:
            y = float(evaluate(x))
            if not np.isfinite(y):
                raise ValueError(f"non-finite signal yield at u2={x}")
            cache[x] = y
        return cache[x]

    peak_index = int(np.nanargmax(signal))
    if 0 < peak_index < len(u2) - 1:
        peak_u2, peak_N = _golden_section_peak(
            value,
            u2[peak_index - 1],
            u2[peak_index + 1],
            log10_tolerance,
        )
    else:
        peak_u2 = float(u2[peak_index])
        peak_N = value(peak_u2)

    if peak_N < N_threshold:
        return {
            "u2_min": np.nan, "u2_max": np.nan,
            "u2_min_open": False, "u2_max_open": False,
            "peak_N": peak_N, "peak_u2": peak_u2,
            "has_sensitivity": False,
        }

    above = signal >= N_threshold
    lower_open = bool(signal[0] >= N_threshold)
    upper_open = bool(signal[-1] >= N_threshold)

    if lower_open:
        lower = np.nan
    else:
        before_peak = np.flatnonzero(above & (u2 <= peak_u2))
        upper_index = int(before_peak[0]) if len(before_peak) else None
        lower_lo = float(u2[upper_index - 1]) if upper_index else float(u2[0])
        lower_hi = float(u2[upper_index]) if upper_index is not None else peak_u2
        if value(lower_hi) < N_threshold:
            lower_hi = peak_u2
        lower = _bisect_log_crossing(
            value,
            lower_lo,
            lower_hi,
            N_threshold,
            rising=True,
            log10_tolerance=log10_tolerance,
        )

    if upper_open:
        upper = np.nan
    else:
        after_peak = np.flatnonzero(above & (u2 >= peak_u2))
        lower_index = int(after_peak[-1]) if len(after_peak) else None
        upper_lo = float(u2[lower_index]) if lower_index is not None else peak_u2
        upper_hi = (
            float(u2[lower_index + 1])
            if lower_index is not None and lower_index + 1 < len(u2)
            else float(u2[-1])
        )
        if value(upper_lo) < N_threshold:
            upper_lo = peak_u2
        upper = _bisect_log_crossing(
            value,
            upper_lo,
            upper_hi,
            N_threshold,
            rising=False,
            log10_tolerance=log10_tolerance,
        )

    return {
        "u2_min": lower, "u2_max": upper,
        "u2_min_open": lower_open, "u2_max_open": upper_open,
        "peak_N": peak_N, "peak_u2": peak_u2,
        "has_sensitivity": True,
    }


def _golden_section_peak(evaluate, lo, hi, log10_tolerance):
    """Maximize a positive-coupling function within a log-space bracket."""
    a, b = np.log10([lo, hi])
    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc = evaluate(10.0 ** c)
    fd = evaluate(10.0 ** d)
    while b - a > log10_tolerance:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = evaluate(10.0 ** c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = evaluate(10.0 ** d)
    x = 10.0 ** ((a + b) / 2.0)
    return x, evaluate(x)


def _bisect_log_crossing(
    evaluate,
    lo,
    hi,
    threshold,
    rising,
    log10_tolerance,
):
    """Bisect a bracketed threshold crossing in log-coupling space."""
    f_lo = evaluate(lo) - threshold
    f_hi = evaluate(hi) - threshold
    bracketed = (f_lo <= 0.0 <= f_hi) if rising else (f_lo >= 0.0 >= f_hi)
    if not bracketed:
        raise ValueError(
            f"unbracketed {'rising' if rising else 'falling'} crossing: "
            f"N({lo})={f_lo + threshold}, N({hi})={f_hi + threshold}"
        )
    a, b = np.log10([lo, hi])
    while b - a > log10_tolerance:
        mid = (a + b) / 2.0
        f_mid = evaluate(10.0 ** mid) - threshold
        if (f_mid >= 0.0) == rising:
            b = mid
        else:
            a = mid
    return 10.0 ** ((a + b) / 2.0)


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
