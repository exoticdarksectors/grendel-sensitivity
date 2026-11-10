"""ATLAS/CMS h -> S S long-lived-particle limits in the BC5 plane.

The PBC 2025 BC5 figure (arXiv:2505.00947, Fig. 22) draws only the limits on
the mixing angle, which BC5 shares with BC4. The LHC searches for Higgs
decays to long-lived particles constrain BR(h -> S S) as a function of
c tau; at the benchmark's fixed BR(h -> S S) = 0.01, the mixing angle only
sets the lifetime (c tau = c tau_0(m_S) / sin^2 theta), so each search
excludes a band of sin^2 theta at every mass where its limit falls below
1 %. The one work that makes that translation is A. Escalante del Valle,
"Dark scalar reinterpretation of searches for Higgs boson decays into
long-lived particles at the LHC", arXiv:2509.02564, JHEP 12 (2025) 129,
Fig. 5; the ESPP 2026 Physics Briefing Book (arXiv:2511.03883, Fig. 8.17)
re-plots its contours. It publishes the bands per search as data:

  github.com/aescalante/LLP_DarkScalarReinterpretation, commit 5788e478,
  plots_sin2theta/combined_exclusion_multi_threshold_sin2theta_log.txt

(one [min, max] excluded sin^2 theta per search and mass point, for
BR(h -> S S) = 1 %, 0.5 % and 0.1 %). This takes the 1 % section as
published: four searches -- CMS muon-system showers (arXiv:2402.01898,
hadronic S decays), CMS displaced dimuons (arXiv:2112.13769 and
arXiv:2402.14491) and the ATLAS/CMS LLP summaries for S -> b b -- each
written as a lower and an upper chain over its mass points, joined by
straight lines in log-log as in the source figure, except across the mass
windows a search vetoes (``GAPS``), where the band is split. The repository
carries no licence; the files are kept for drawing with the citation.

Output: ``curves/data/bc5/existing/lhc_hss_*.dat`` (``mass_GeV sin2theta``,
for each piece of the band the lower chain, a blank line, the upper chain).

    python -m curves.tools.lhc_bc5 [--source FILE]
"""
from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from pathlib import Path

from ..plot.data import BC5_DATA

DATE = "2026-09-29"
REPO = "aescalante/LLP_DarkScalarReinterpretation"
COMMIT = "5788e478a637c692ba5500c6aed90959ec8fb4dc"
FILE = "plots_sin2theta/combined_exclusion_multi_threshold_sin2theta_log.txt"
URL = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{FILE}"
OUT_DIR = BC5_DATA / "existing"
SCENARIO = r"= 1\%"
SEARCHES = {   # source label fragment -> (output stem, description)
    "LHC LLP Summaries": ("lhc_hss_bb_llp_summaries", "S -> b b, ATLAS/CMS LLP summary plots (ATL-PHYS-PUB-2025-002)"),
    "arXiv:2402.14491": ("lhc_hss_mumu_2402.14491", "S -> mu mu, CMS displaced dimuons, arXiv:2402.14491"),
    "arXiv:2112.13769": ("lhc_hss_mumu_2112.13769", "S -> mu mu, CMS displaced dimuons, arXiv:2112.13769"),
    "arXiv:2402.01898": ("lhc_hss_hadronic_2402.01898",
                         "S -> pi pi, K K, g g, s s, c c, tau_h tau_h, CMS muon-system showers, arXiv:2402.01898"),
}
# Mass windows a search vetoes: its points stop at the window edges and the
# band must not be drawn across them. arXiv:2112.13769 vetoes the J/psi and
# psi(2S) (2.91-3.89 GeV) and the Upsilon(1S-3S) (8.99-10.78 GeV) dimuon
# windows; the ESPP 2026 Briefing Book (Fig. 8.17) draws the band with these
# gaps.
GAPS = {"arXiv:2112.13769": ((2.91, 3.89), (8.99, 10.78))}


def pieces(rows, gaps):
    """``rows`` split wherever two consecutive mass points span a vetoed window."""
    out = [[rows[0]]]
    for prev, row in zip(rows, rows[1:]):
        if any(abs(prev[0] - a) < 1e-6 and abs(row[0] - b) < 1e-6 for a, b in gaps):
            out.append([])
        out[-1].append(row)
    return out


def parse(text: str) -> dict[str, list[tuple[float, float, float]]]:
    """``{label: [(m, min, max), ...]}`` of the BR(h -> S S) = 1 % section."""
    sections = re.split(r"(?m)^# Scenario: ", text)
    wanted = [s for s in sections if s.splitlines() and SCENARIO in s.splitlines()[0]]
    if len(wanted) != 1:
        raise RuntimeError(f"expected one BR(h->SS) = 1% section, found {len(wanted)}")
    out: dict[str, list[tuple[float, float, float]]] = {}
    label = None
    for line in wanted[0].splitlines()[1:]:
        line = line.strip()
        if line.startswith("# Label:"):
            label = line[len("# Label:"):].strip()
            out[label] = []
        elif line and not line.startswith(("#", "m(S)")) and label is not None:
            m, lo, hi = (float(v) for v in line.split(","))
            out[label].append((m, lo, hi))
    return out


def write(bands: dict[str, list[tuple[float, float, float]]]) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    for fragment, (stem, description) in SEARCHES.items():
        matches = [label for label in bands if fragment in label]
        if len(matches) != 1:
            raise RuntimeError(f"{fragment}: expected one search, found {matches}")
        rows = sorted(bands[matches[0]])
        if any(not (0.0 < lo < hi) for _, lo, hi in rows):
            raise RuntimeError(f"{fragment}: an empty or inverted band")
        path = OUT_DIR / f"{stem}.dat"
        with path.open("w") as fh:
            fh.write(f"# LHC h -> S S limit in the BC5 plane at BR(h -> S S) = 0.01: {description}\n")
            fh.write("# Source: A. Escalante del Valle, arXiv:2509.02564, JHEP 12 (2025) 129, Fig. 5;\n")
            fh.write(f"#   data github.com/{REPO} at {COMMIT[:8]}, {FILE}\n")
            fh.write(f"#   (label {matches[0]!r}), taken as published.\n")
            gaps = GAPS.get(fragment, ())
            fh.write("# The excluded band lies between each pair of chains (lower, then upper), over the\n")
            fh.write("#   source's mass points, joined by straight lines in log-log.\n")
            if gaps:
                windows = ", ".join(f"{a:g}-{b:g}" for a, b in gaps)
                fh.write(f"# Split at the vetoed mass windows {windows} GeV (no exclusion there).\n")
            fh.write(f"# Extracted: {DATE}\n")
            fh.write("# mass_GeV  sin2theta\n")
            for k, piece in enumerate(pieces(rows, gaps)):
                if k:
                    fh.write("\n")
                for m, lo, _ in piece:
                    fh.write(f"{m:.6g}  {lo:.10e}\n")
                fh.write("\n")
                for m, _, hi in piece:
                    fh.write(f"{m:.6g}  {hi:.10e}\n")
        written.append(path)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=None, help="a local copy of the data file (default: fetch it)")
    args = ap.parse_args(argv)
    if args.source is not None:
        text = args.source.read_text()
    else:
        with urllib.request.urlopen(URL, timeout=60) as response:
            text = response.read().decode()
    for path in write(parse(text)):
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
