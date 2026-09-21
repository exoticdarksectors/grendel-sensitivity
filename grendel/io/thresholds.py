"""Project a secondary threshold out of a sensitivity table.

    python -m grendel.io.thresholds RESULTS/bc4/sensitivity.csv --threshold 10

A scan run with ``--thresholds 3 10`` solves every yield curve twice and
keeps the second island in ``<field>_N10`` columns next to the primary
ones. The renderer reads the second island from its own table,
``sensitivity_nsig10.csv``, with the ordinary column names; this writes
that table by moving the ``_N<T>`` columns into place and dropping the
primary island. The two tables then describe the identical Monte Carlo
realisation solved at two thresholds, which is what lets the dashed
contour be drawn inside the filled one.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def threshold_tag(threshold: float) -> str:
    return f"N{threshold:g}"


def split_threshold(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """The rows of ``frame`` with the ``_N<threshold>`` island in place of the
    primary one. Columns without a secondary counterpart are kept as they are."""
    suffix = "_" + threshold_tag(threshold)
    secondary = [c for c in frame.columns if c.endswith(suffix)]
    if not secondary:
        raise ValueError(f"no {suffix} columns: the scan was not solved at N >= {threshold:g}")
    primary = [c[: -len(suffix)] for c in secondary]
    other_suffixes = {c[c.rindex("_N"):] for c in frame.columns
                      if "_N" in c and c[c.rindex("_N") + 2:].replace(".", "").isdigit()} - {suffix}
    out = frame.drop(columns=primary, errors="ignore")
    out = out.drop(columns=[c for c in out.columns if any(c.endswith(s) for s in other_suffixes)])
    out = out.rename(columns=dict(zip(secondary, primary)))
    # Keep the primary column order where it exists.
    ordered = [c for c in frame.columns if c in out.columns]
    return out[ordered]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sensitivity", type=Path, help="sensitivity.csv of a scan run with several thresholds")
    ap.add_argument("--threshold", type=float, default=10.0)
    ap.add_argument("--out", type=Path, default=None,
                    help="output table (default: sensitivity_nsig<T>.csv beside the input)")
    args = ap.parse_args(argv)
    frame = pd.read_csv(args.sensitivity)
    out = args.out or args.sensitivity.with_name(f"sensitivity_nsig{args.threshold:g}.csv")
    split_threshold(frame, args.threshold).to_csv(out, index=False)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
