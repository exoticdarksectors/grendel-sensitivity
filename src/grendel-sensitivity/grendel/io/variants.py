"""Project one variant curve out of a sensitivity table.

A scan may solve extra exclusion islands on the same Monte Carlo besides the nominal one -- BC5 solves
every production mode alone and the curves rescaled to other values of BR(h -> SS) -- and writes them
as ``<variant>_<field>`` columns next to the nominal ``<field>`` ones (``hSS_u2_min``, ``hSS_peak_N``,
``brhss0.001_u2_max_N10``, ...). This module rewrites such a table with one variant in the nominal
place, in the layout of ``sensitivity.csv``, so that whatever plots ``sensitivity.csv`` plots the
variant.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ANCHOR_FIELDS = ("u2_min", "invf_min")


def variants_in(frame: pd.DataFrame) -> list[str]:
    """The variant prefixes of ``frame``: every ``X`` with an ``X_<anchor>`` column, where the anchor
    is the table's lower island edge (``u2_min`` or ``invf_min``)."""
    for anchor in ANCHOR_FIELDS:
        if anchor in frame.columns:
            suffix = "_" + anchor
            return [c[: -len(suffix)] for c in frame.columns if c.endswith(suffix)]
    return []


def split_variant(frame: pd.DataFrame, variant: str) -> pd.DataFrame:
    """The rows of ``frame`` with the ``<variant>_*`` island in place of the nominal one; the nominal
    columns it replaces and every other variant's copies of them are dropped."""
    head = variant + "_"
    renamed = [c for c in frame.columns if c.startswith(head)]
    if not renamed:
        known = ", ".join(variants_in(frame)) or "none"
        raise ValueError(f"no {head}* columns in the table (variants present: {known})")
    bases = [c[len(head):] for c in renamed]
    base_set = set(bases)
    dropped = [c for c in frame.columns
               if c not in renamed and (c in base_set or any(c.endswith("_" + b) for b in bases))]
    out = frame.drop(columns=dropped).rename(columns=dict(zip(renamed, bases)))
    order = []
    for c in frame.columns:
        if c in renamed:
            continue
        if c in dropped:
            if c in base_set:
                order.append(c)
            continue
        order.append(c)
    order.extend(b for b in bases if b not in order)
    return out[order]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sensitivity", type=Path, help="sensitivity.csv of a scan that solved variant curves")
    ap.add_argument("--variant", default=None,
                    help="the variant to project out (its column prefix, e.g. hSS, mixing, brhss0.001)")
    ap.add_argument("--list", action="store_true", help="print the variants of the table and exit")
    ap.add_argument("--out", type=Path, default=None,
                    help="output table (default: sensitivity_<variant>.csv beside the input)")
    args = ap.parse_args(argv)
    frame = pd.read_csv(args.sensitivity)
    if args.list or args.variant is None:
        for name in variants_in(frame):
            print(name)
        if args.variant is None and not args.list:
            ap.error("--variant is required (--list shows the variants of the table)")
        return 0
    out = args.out or args.sensitivity.with_name(f"{args.sensitivity.stem}_{args.variant}.csv")
    split_variant(frame, args.variant).to_csv(out, index=False)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
