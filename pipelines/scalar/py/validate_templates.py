"""Compare regenerated exHad scalar-1809 bundles with the published ones (statistically: Linux vs macOS arm64).

    python validate_templates.py OURS_DIR REF_DIR REPORT.tsv

Exact fields: mass, ctau_m_u2eq1 (exHad's table), n_templates, base seed, per-point exHad seed. Statistical:
channel-group fractions (two-sample binomial z), mean multiplicity and charged multiplicity (z of the difference
of means), the fraction with >= 2 charged daughters (the visibility the selection needs), and two-sample KS tests
of the daughter |p| and of the charged-daughter |p|. Byte identity of the physics arrays is reported too.
Exit status 1 if an exact field differs or a statistical test fails at |z| > 5 or KS p < 1e-6 (global).
"""
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ours, ref, report = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
REPRO = Path(__file__).resolve().parents[3]  # the report path is printed relative to the repository
PHYS = ("daughter_counts", "pdg", "px", "py", "pz", "energy", "mass", "charge", "channel_label", "n_charged")


def z_binom(a, na, b, nb):
    p = (a + b) / (na + nb)
    s = np.sqrt(max(p * (1 - p), 1e-300) * (1 / na + 1 / nb))
    return 0.0 if s == 0 else (a / na - b / nb) / s


def z_mean(x, y):
    s = np.sqrt(x.var(ddof=1) / len(x) + y.var(ddof=1) / len(y))
    return 0.0 if s == 0 else (x.mean() - y.mean()) / s


rows, bad = [], []
files = sorted(ref.glob("templates_*.npz"))
for f in files:
    g = ours / f.name
    if not g.exists():
        bad.append(f"{f.name}: missing")
        continue
    a, b = np.load(g), np.load(f)
    m = float(b["mass_GeV"])
    exact = {k: (a[k].item() == b[k].item()) for k in ("mass_GeV", "ctau_m_u2eq1", "n_templates", "seed", "exhad_seed")}
    identical = all(np.array_equal(a[k], b[k]) for k in PHYS)
    la, lb = a["channel_label"].astype(str), b["channel_label"].astype(str)
    na, nb = len(la), len(lb)
    zch = {c: z_binom((la == c).sum(), na, (lb == c).sum(), nb) for c in sorted(set(la) | set(lb))}
    zmult = z_mean(a["daughter_counts"].astype(float), b["daughter_counts"].astype(float))
    zchg = z_mean(a["n_charged"].astype(float), b["n_charged"].astype(float))
    zvis = z_binom((a["n_charged"] >= 2).sum(), na, (b["n_charged"] >= 2).sum(), nb)
    # |p| to 9 significant digits: two-body decays are monochromatic, and last-bit jitter would make the KS
    # statistic of two delta functions meaningless
    pa = np.array(
        [float(f"{x:.9g}") for x in np.sqrt(a["px"] ** 2 + a["py"] ** 2 + a["pz"] ** 2)])
    pb = np.array([float(f"{x:.9g}") for x in np.sqrt(b["px"] ** 2 + b["py"] ** 2 + b["pz"] ** 2)])
    # event by event: same daughters (pdg) and momenta equal to 1e-6 relative (exHad seeds every sub-request, so
    # most events are the same decays up to floating-point noise)
    same_ev = np.nan
    if np.array_equal(a["daughter_counts"], b["daughter_counts"]):
        off = np.concatenate([[0], np.cumsum(b["daughter_counts"])])
        e = b["energy"]
        close = (a["pdg"] == b["pdg"]) & np.all(
            [np.abs(a[k] - b[k]) <= 1e-6 * np.maximum(e, 1e-9) for k in ("px", "py", "pz", "energy")], axis=0)
        same_ev = float(np.logical_and.reduceat(close, off[:-1]).mean())
    ks = stats.ks_2samp(pa, pb).pvalue if not identical else 1.0
    ca, cb = a["charge"] != 0, b["charge"] != 0
    ksc = stats.ks_2samp(pa[ca], pb[cb]).pvalue if (ca.any() and cb.any() and not identical) else 1.0
    zmax = max([abs(v) for v in zch.values()] + [abs(zmult), abs(zchg), abs(zvis)])
    rows.append((m, f.name, all(exact.values()), identical, same_ev, zmax, ks, ksc,
                 float(a["daughter_counts"].mean()), float(b["daughter_counts"].mean()),
                 float((a["n_charged"] >= 2).mean()), float((b["n_charged"] >= 2).mean()),
                 ";".join(f"{c}:{(la == c).mean():.4f}/{(lb == c).mean():.4f}" for c in zch)))
    if not all(exact.values()):
        bad.append(f"{f.name}: exact fields differ {[k for k, v in exact.items() if not v]}")
    if zmax > 5:
        bad.append(f"{f.name}: |z| = {zmax:.1f}")
    if min(ks, ksc) < 1e-6:
        bad.append(f"{f.name}: KS p = {min(ks, ksc):.2g}")

hdr = ("mass_GeV\tfile\texact_fields_equal\tphysics_byte_identical\tfrac_events_equal_1e-6\tmax_abs_z\tKS_p_all\tKS_p_charged\t"
       "mult_ours\tmult_ref\tvis2_ours\tvis2_ref\tchannel_fractions_ours/ref")
report.write_text(hdr + "\n" + "\n".join("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r)
                                         for r in rows) + "\n")
z = np.array([r[5] for r in rows])
ks = np.array([min(r[6], r[7]) for r in rows])
se = np.array([r[4] for r in rows], float)
print(f"{len(rows)}/{len(files)} bundles compared; exact fields equal in {sum(r[2] for r in rows)}; "
      f"physics arrays byte-identical in {sum(r[3] for r in rows)}")
print(f"event by event (same daughters, momenta within 1e-6): fraction of events equal per bundle: median "
      f"{np.nanmedian(se):.4f}, min {np.nanmin(se):.4f}; daughter-count arrays differ in {int(np.isnan(se).sum())} bundles")
print(f"max |z| over channel fractions, <mult>, <n_charged>, P(>=2 charged): median {np.median(z):.2f}, "
      f"max {z.max():.2f}; KS p-values: min {ks.min():.3g}, fraction < 0.01: {(ks < 0.01).mean():.3f}")
print(f"report: {report.resolve().relative_to(REPRO) if report.resolve().is_relative_to(REPRO) else report}")
for line in bad:
    print("FAIL", line)
sys.exit(1 if bad else 0)
