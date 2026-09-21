"""Primitives of the single-source variation envelope (BC4, BC10).

Every uncertainty campaign varies one source at a time -- FONLL scale
choices, NNPDF replicas, heavy-quark masses, decay-model alternates -- and
re-runs the whole chain for each. The combination works in
``x = log10(coupling)`` because an exclusion edge spans decades:

* a *named* axis (scale, mass, decay model) contributes the interval between
  the extreme members and the central, :func:`axis_interval`;
* the PDF replica ensemble contributes its 16th and 84th percentiles,
  :func:`pdf_percentiles` (the raw extrema of 100 replicas are deliberately
  not a headline interval);
* the display envelope is the outermost endpoint over the per-source
  intervals. It is neither a quadrature combination nor a confidence band.

Same-physics repeats with fresh seeds are *numerical controls*: their
spread, :func:`repeat_statistics`, is reported next to the physics intervals
and never enters the envelope.

The campaign's central curve carries campaign statistics; the published
central curve may be a higher-statistics run. A band is therefore *rebased*
onto the canonical central: in dex (:func:`rebase_dex`, BC4) or
multiplicatively, which is the same shift written on the linear value
(BC10). Where the two models genuinely differ -- which intervals include the
central endpoint, which columns are rebased, how a missing member is
reported -- lives in their own ``combine`` modules, next to their column
contracts.

The HNL production band uses a different rule (:mod:`.quadrature`).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def log10_boundary(value) -> float | None:
    """``log10(value)`` for a finite, positive boundary; None otherwise."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value) or value <= 0:
        return None
    return float(np.log10(value))


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float
    lo_name: str
    hi_name: str

    def magnitude(self, xc: float) -> float:
        """The larger displacement of either endpoint from ``xc``."""
        return max(abs(self.lo - xc), abs(self.hi - xc))


def axis_interval(members, *, central: tuple[str, float] | None = None) -> Interval:
    """The extrema over ``members`` ``[(name, x), ...]``, optionally including
    the central point (BC4 and BC10 both span every named axis together with
    the central)."""
    candidates = list(members)
    if central is not None:
        candidates = [central, *candidates]
    if not candidates:
        raise ValueError("axis_interval needs at least one member")
    lo_name, lo = min(candidates, key=lambda item: item[1])
    hi_name, hi = max(candidates, key=lambda item: item[1])
    return Interval(float(lo), float(hi), lo_name, hi_name)


@dataclass(frozen=True)
class PdfStatistics:
    p16: float
    p84: float
    std: float
    n: int


def pdf_percentiles(xs) -> PdfStatistics:
    """16th/84th percentiles (linear interpolation) and sample standard
    deviation (``ddof=1``, NaN below two members) of replica boundaries.
    With no members the percentiles are NaN; with one they are that value."""
    xs = np.asarray(list(xs), dtype=float)
    if len(xs) >= 2:
        p16, p84 = np.quantile(xs, [0.16, 0.84], method="linear")
        return PdfStatistics(float(p16), float(p84), float(np.std(xs, ddof=1)), len(xs))
    if len(xs) == 1:
        return PdfStatistics(float(xs[0]), float(xs[0]), np.nan, 1)
    return PdfStatistics(np.nan, np.nan, np.nan, 0)


@dataclass(frozen=True)
class RepeatStatistics:
    median_abs_dex: float
    max_abs_dex: float
    median_abs_fraction: float
    max_abs_fraction: float
    n: int


def repeat_statistics(xc: float, repeats) -> RepeatStatistics:
    """Absolute boundary shifts of same-physics repeats relative to the
    campaign central, in dex and as a linear fraction."""
    xs = np.asarray(list(repeats), dtype=float)
    if len(xs) == 0:
        return RepeatStatistics(np.nan, np.nan, np.nan, np.nan, 0)
    abs_dex = np.abs(xs - xc)
    abs_fraction = np.abs(10.0 ** (xs - xc) - 1.0)
    return RepeatStatistics(float(np.median(abs_dex)), float(np.max(abs_dex)),
                            float(np.median(abs_fraction)), float(np.max(abs_fraction)),
                            len(xs))


def rebase_dex(xref: float, x: float, xc: float) -> float:
    """Apply the shift ``x - xc`` measured in the campaign to the canonical
    central ``xref``, back on the linear scale."""
    return 10.0 ** (xref + x - xc)
