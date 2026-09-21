"""Per-experiment projections from each experiment's own publication.

The compilation figures (the PBC BC7 panel, HNLimits) redraw contours
that first appeared in the experiments' papers. This command records, for
ANUBIS, MATHUSLA, CODEX-b and SHiP, what the original source is and what
"original" can mean for each:

* SHiP (arXiv:1811.00930) and MATHUSLA (arXiv:1806.07396, digitised once
  by Bondarenko et al. arXiv:1805.08567) have no standalone data release;
  the HNLimits-bundled files ARE the original-source data, and are copied
  alongside their provenance note.
* CODEX-b (arXiv:1911.00481) publishes no numerical contour; the curve is
  vector-extracted from the EoI's own per-flavour panels by
  ``vector_codexb_hnl``.
* ANUBIS (arXiv:2606.26862 official EW production; arXiv:2512.13011
  Wang-Zhang meson production) likewise: the tracked ``.dat`` files are
  extractions from those figures.

It caches the arXiv PDFs under ``<work dir>/curves/original_sources_pdfs``
and writes ``SOURCE.md`` per experiment under
``curves/data/hnl/original_sources/``.

    python -m curves.tools.original_sources
"""
from __future__ import annotations

import os
import shutil
import sys
import urllib.request
from pathlib import Path

from ..plot import hnl_panel
from ..plot.data import HNL_DATA

OUT = HNL_DATA / "original_sources"
PDF_DIR = Path(os.environ.get("GRENDEL_WORK_DIR", "grendel_work")) / "curves" / "original_sources_pdfs"

SOURCES = [
    {"exp": "SHiP", "arxiv": ["1811.00930"], "figure": "Fig. 33 in 1811.00930",
     "note": "Original SHiP HNL projection. Data is identical to the HNLimits-bundled copy under SHiP/ and "
             "1811.00930/; copied here.",
     "hnl_files": ["SHiP/Umu4.dat", "SHiP/Umu4_top.dat", "1811.00930/Umu4_optimistic.dat",
                   "1811.00930/Umu4_conservative.dat"]},
    {"exp": "MATHUSLA", "arxiv": ["1806.07396", "1903.04497", "1805.08567"],
     "figure": "Fig. 7 in 1805.08567 (digitized once from MATHUSLA paper)",
     "note": "Bondarenko et al. is the de-facto canonical copy. MATHUSLA40 (updated geometry) numbers are not "
             "yet public.",
     "hnl_files": ["1805.08567/mathusla_lower.dat", "1805.08567/mathusla_upper.dat"]},
    {"exp": "CODEX-b", "arxiv": ["1911.00481", "2203.07316"],
     "figure": "figures/HNL_{Ele,Mu,Tau}.pdf of 1911.00481v2 (fig:{Ne,mu,tau}_combined)",
     "note": "No public machine-readable data. curves.tools.vector_codexb_hnl vector-extracts the CODEX-b "
             "contour from the EoI arXiv source's own three single-flavour panels -- all of Ue/Umu/Utau, "
             "unlike the muon-only PBC BC7 redraw -- and converts the EoI Dirac convention to Majorana. "
             "The 2025 ESPP update arXiv:2505.05952 carries no HNL contour of its own and defers to this EoI.",
     "hnl_files": [], "extracted": "curves.tools.vector_codexb_hnl",
     "extracted_files": ["CODEX-b_2019_Ue.dat", "CODEX-b_2019_Umu.dat", "CODEX-b_2019_Utau.dat"]},
    {"exp": "ANUBIS", "arxiv": ["2001.04750", "2512.13011", "2606.26862"],
     "figure": "Fig. 4 of 2606.26862 for official EW production; Fig. 1 of 2512.13011 for Wang-Zhang meson "
               "production",
     "note": "Official 2026 ANUBIS HNL contours are EW-production driven for Ue/Umu. Wang-Zhang 2512.13011 is "
             "electron-only and charm/bottom meson-production driven; muon is only stated to be similar apart "
             "from threshold effects.",
     "hnl_files": [], "extracted": "vector extraction from the figures above",
     "extracted_files": ["ANUBIS_2026_Ue.dat", "ANUBIS_2026_Umu.dat", "ANUBIS_WangZhang_2025_ceiling_Ue.dat"]},
]


def fetch_arxiv(arxiv_id, dest_dir):
    out = dest_dir / f"{arxiv_id}.pdf"
    if out.exists():
        return out, "cached"
    url = f"https://arxiv.org/pdf/{arxiv_id}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "grendel-sensitivity/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r, out.open("wb") as fh:
            fh.write(r.read())
        return out, "downloaded"
    except Exception as e:  # noqa: BLE001 - a failed download is reported, not fatal
        return None, f"failed: {e}"


def main() -> int:
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    hnl_data = hnl_panel.hnlimits_data_root()
    for src in SOURCES:
        exp_dir = OUT / src["exp"]
        exp_dir.mkdir(parents=True, exist_ok=True)
        status = {aid: fetch_arxiv(aid, PDF_DIR) for aid in src["arxiv"]}
        copied = []
        for rel in src["hnl_files"]:
            srcfile = hnl_data / rel
            if srcfile.exists():
                dst = exp_dir / Path(rel).name
                shutil.copyfile(srcfile, dst)
                copied.append(dst.relative_to(HNL_DATA))
        with (exp_dir / "SOURCE.md").open("w") as fh:
            fh.write(f"# {src['exp']} — original-source provenance\n\n")
            fh.write(f"**Citations:** {', '.join('arXiv:' + a for a in src['arxiv'])}\n\n")
            fh.write(f"**Figure containing the contour:** {src['figure']}\n\n")
            fh.write(f"**Note:** {src['note']}\n\n")
            fh.write("**arXiv PDFs cached at:**\n")
            for aid, (path, state) in status.items():
                fh.write(f"- {aid}: {state}" + (f" -> {path}" if path else "") + "\n")
            if copied:
                fh.write("\n**Copied from the HNLimits bundle (these *are* the original-source data):**\n")
                for c in copied:
                    fh.write(f"- curves/data/hnl/{c}\n")
            elif src.get("extracted"):
                fh.write(f"\n**Data produced by `{src['extracted']}`:**\n")
                for name in src.get("extracted_files", []):
                    fh.write(f"- curves/data/hnl/original_sources/{src['exp']}/{name}\n")
        print(f"{src['exp']:10s} {', '.join(f'{a}: {s[1]}' for a, s in status.items())},  copied: {len(copied)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
