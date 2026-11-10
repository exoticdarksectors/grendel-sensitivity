"""Compare the exHad ALP-fermion templates of step 03 with the stored macOS bundles of bc10_v24_20261019.

exHad's same-seed event streams are not bit-identical between Linux/x86-64/gcc and macOS/arm64/clang (libm last-bit
differences; see the shared notes), so per mass the check is statistical, on two samples of 20 000 decays each:

  exact      ctau_m_u2eq1, ctau_m_u2eq1_reference, ctau_convention_factor, ctau_m_u2eq1_exhad_table, exhad_seed,
             n_templates (table values, must be bit-identical);
  event      the two streams are compared decay by decay: whether the PDG sequence is the same, the fraction of
             decays bit-identical in (pdg, px, py, pz, E), and the largest |delta p_i| or |delta E| of any daughter [GeV];
  chi2       final-state classes (the sorted PDG list of a decay; the 12 most common + rest), exHad channel labels,
             number of charged daughters, number of daughters: Pearson chi2 test of homogeneity, p-value;
  KS         |p| of the daughters per species (pi+-, K+-, mu, e, gamma, K_L, p/n), two-sample KS p-value, on |p|
             rounded to 1e-6 GeV (two-body decays at rest give one |p| value, which last-bit noise would otherwise split);
  numbers    visible fraction (>= 2 charged), mean multiplicity, mean charged |p|: difference / its MC error.

A mass passes when every exact field matches, and is flagged when a p-value is below 1e-3 or a pull exceeds 4
(with ~15 tests per mass and 151 masses, a few p < 0.01 are expected by chance). When the PDG sequences agree the
templates are the same decays and the momentum difference is the whole story.

    python validate_templates.py <our template dir> <report.txt> [mass ...]

The stored bundles are reference material ($GRENDEL_REF/grendel/campaigns/bc10/bc10_v24_20261019/templates_exhad);
without them the script prints SKIP and exits 0.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from ref import REF as REF_TREE, skip_unless  # noqa: E402  (pipelines/common/ref.py)

REF = REF_TREE / "grendel/campaigns/bc10/bc10_v24_20261019/templates_exhad"   # $GRENDEL_REF (optional)
EXACT = ("ctau_m_u2eq1", "ctau_m_u2eq1_reference", "ctau_convention_factor", "ctau_m_u2eq1_exhad_table",
         "exhad_seed", "n_templates", "includes_full_branching", "decay_backend")
SPECIES = {"pi+-": (211,), "K+-": (321,), "mu": (13,), "e": (11,), "gamma": (22,), "K_L": (130,),
           "p/n": (2212, 2112)}


def load(path):
    z = np.load(path)
    d = {k: z[k] for k in z.files}
    off = np.concatenate([[0], np.cumsum(d["daughter_counts"])])
    d["offsets"] = off
    d["p"] = np.sqrt(d["px"] ** 2 + d["py"] ** 2 + d["pz"] ** 2)
    d["p_rounded"] = np.round(d["p"], 6)
    pdg = d["pdg"]
    d["classes"] = [tuple(sorted(pdg[off[i]:off[i + 1]].tolist())) for i in range(len(off) - 1)]
    return d


def chi2_p(a, b, top=None):
    ca, cb = Counter(a), Counter(b)
    keys = [k for k, _ in (ca + cb).most_common(top)]
    rows = [[ca.get(k, 0) for k in keys], [cb.get(k, 0) for k in keys]]
    if top is not None:
        rows[0].append(sum(ca.values()) - sum(rows[0]))
        rows[1].append(sum(cb.values()) - sum(rows[1]))
    t = np.array(rows, dtype=float)
    t = t[:, t.sum(axis=0) > 0]
    if t.shape[1] < 2:
        return 1.0
    return float(stats.chi2_contingency(t, correction=False)[1])


def identical_fraction(a, b):
    n = min(len(a["daughter_counts"]), len(b["daughter_counts"]))
    same = 0
    for i in range(n):
        sa, ea = a["offsets"][i], a["offsets"][i + 1]
        sb, eb = b["offsets"][i], b["offsets"][i + 1]
        if ea - sa != eb - sb:
            continue
        if (np.array_equal(a["pdg"][sa:ea], b["pdg"][sb:eb])
                and all(np.array_equal(a[k][sa:ea], b[k][sb:eb]) for k in ("px", "py", "pz", "energy"))):
            same += 1
    return same / n


def pull(x, y, sx, sy):
    """(x - y) / error; 0 for a difference below 1e-9 relative (floating-point noise of identical decays, e.g. a
    two-body decay at rest, whose |p| has no spread)."""
    if abs(x - y) <= 1e-9 * max(abs(x), abs(y)):
        return 0.0
    s = np.hypot(sx, sy)
    return float((x - y) / s) if s > 0 else float("inf")


def compare(ours, ref):
    a, b = load(ours), load(ref)
    res = {"exact_bad": [k for k in EXACT if k in b and not np.array_equal(a.get(k), b[k])]}
    res["identical"] = identical_fraction(a, b)
    res["same_pdg"] = (np.array_equal(a["daughter_counts"], b["daughter_counts"])
                       and np.array_equal(a["pdg"], b["pdg"]))
    res["max_dp"] = (float(np.max(np.abs(np.stack([a[k] - b[k] for k in ("px", "py", "pz", "energy")]))))
                     if res["same_pdg"] else float("nan"))
    pv = {"class": chi2_p(a["classes"], b["classes"], top=12),
          "label": chi2_p(a["channel_label"].tolist(), b["channel_label"].tolist()),
          "n_charged": chi2_p(a["n_charged"].tolist(), b["n_charged"].tolist()),
          "n_daughters": chi2_p(a["daughter_counts"].tolist(), b["daughter_counts"].tolist())}
    for name, codes in SPECIES.items():
        sa = a["p_rounded"][np.isin(np.abs(a["pdg"]), codes)]
        sb = b["p_rounded"][np.isin(np.abs(b["pdg"]), codes)]
        if len(sa) >= 20 and len(sb) >= 20:
            pv[f"KS p({name})"] = float(stats.ks_2samp(sa, sb).pvalue)
    res["p"] = pv
    n = len(a["daughter_counts"])
    va, vb = (a["n_charged"] >= 2).mean(), (b["n_charged"] >= 2).mean()
    ma, mb = a["daughter_counts"].mean(), b["daughter_counts"].mean()
    qa, qb = a["p"][a["charge"] != 0], b["p"][b["charge"] != 0]
    res["vis"] = (va, vb, pull(va, vb, np.sqrt(va * (1 - va) / n), np.sqrt(vb * (1 - vb) / n)))
    res["mult"] = (ma, mb, pull(ma, mb, a["daughter_counts"].std() / np.sqrt(n), b["daughter_counts"].std() / np.sqrt(n)))
    res["pch"] = (qa.mean(), qb.mean(), pull(qa.mean(), qb.mean(), qa.std() / np.sqrt(len(qa)), qb.std() / np.sqrt(len(qb))))
    res["ctau"] = (float(a["ctau_m_u2eq1"]), float(b["ctau_m_u2eq1"]))
    return res


def main(our_dir, report, *masses):
    skip_unless("BC10 templates vs the stored macOS bundles", REF)
    our_dir = Path(our_dir)
    files = sorted(our_dir.glob("templates_*.npz"))
    if masses:
        want = {f"templates_{float(m):.3f}".replace(".", "p") + ".npz" for m in masses}
        files = [f for f in files if f.name in want]
    lines, n_exact_bad, flagged, allp, n_same, max_dp = [], 0, [], [], 0, []
    for f in files:
        ref = REF / f.name
        if not ref.exists():
            lines.append(f"{f.name}: no stored bundle")
            continue
        r = compare(f, ref)
        allp += list(r["p"].values())
        n_same += r["same_pdg"]
        max_dp.append(r["max_dp"])
        n_exact_bad += bool(r["exact_bad"])
        minp = min(r["p"].values())
        pulls = [abs(r[k][2]) for k in ("vis", "mult", "pch")]
        flag = bool(r["exact_bad"]) or minp < 1e-3 or max(pulls) > 4
        if flag:
            flagged.append(f.name)
        lines.append(
            f"{f.stem[10:]}: exact {'OK' if not r['exact_bad'] else 'DIFFER ' + ','.join(r['exact_bad'])}"
            f" (ctau {r['ctau'][0]:.10e})  same PDG sequence {r['same_pdg']}, max |dp| {r['max_dp']:.1e} GeV, "
            f"bit-identical decays {r['identical']:.3f}  "
            f"vis {r['vis'][0]:.4f}/{r['vis'][1]:.4f} ({r['vis'][2]:+.1f})  "
            f"mult {r['mult'][0]:.3f}/{r['mult'][1]:.3f} ({r['mult'][2]:+.1f})  "
            f"<p_ch> {r['pch'][0]:.4f}/{r['pch'][1]:.4f} ({r['pch'][2]:+.1f})  min p {minp:.3f}"
            + ("  FLAG" if flag else "") + "\n    " + "  ".join(f"{k} {v:.3f}" for k, v in r["p"].items()))
    if not files:
        print("no template bundles to compare")
        return 1
    allp = np.array(allp)
    head = (f"templates compared: {len(files)}; table fields bit-identical in {len(files) - n_exact_bad}; "
            f"same decays (PDG sequence) in {n_same}, largest daughter |dp| or |dE| {np.nanmax(max_dp) if np.isfinite(max_dp).any() else float('nan'):.1e} GeV; "
            f"{len(allp)} homogeneity tests: fraction p<0.05 {np.mean(allp < 0.05):.3f}, p<0.01 {np.mean(allp < 0.01):.3f}, "
            f"min p {allp.min() if len(allp) else float('nan'):.2e} (uniform expectation 0.05 / 0.01); "
            f"flagged: {len(flagged)} {' '.join(flagged)}")
    Path(report).write_text(head + "\n\n" + "\n".join(lines) + "\n")
    print(head)
    return 0 if n_exact_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
