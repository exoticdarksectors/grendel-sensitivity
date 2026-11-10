"""Compare two exHad template bundles statistically (independent 20k-decay samples of the same mass).

    python compare_templates.py <new.npz> <published.npz>

Prints the lifetime, the chi2 of the channel-label and charged-multiplicity distributions, and two-sample KS
p-values of the charged-daughter |p| and of the visible (charged) momentum sum per decay.
"""
import sys

import numpy as np
from scipy.stats import chi2 as C2, ks_2samp


def load(p):
    d = np.load(p)
    n = d["daughter_counts"]
    owner = np.repeat(np.arange(len(n)), n)
    st = d["stable"].astype(bool)
    ch = st & (d["charge"] != 0)
    p3 = np.sqrt(d["px"] ** 2 + d["py"] ** 2 + d["pz"] ** 2)
    vis = np.zeros((len(n), 3))
    for k, c in enumerate(("px", "py", "pz")):
        np.add.at(vis[:, k], owner[ch], d[c][ch])
    return dict(ctau=float(d["ctau_m_u2eq1"]), labels=d["channel_label"], nch=d["n_charged"],
                p_ch=p3[ch], pvis=np.linalg.norm(vis, axis=1), mult=n)


def chi2_hist(a, b):
    keys = sorted(set(a) | set(b))
    ca = np.array([np.sum(a == k) for k in keys], float); cb = np.array([np.sum(b == k) for k in keys], float)
    na, nb = ca.sum(), cb.sum()
    use = (ca + cb) > 0
    x = (ca[use] * np.sqrt(nb / na) - cb[use] * np.sqrt(na / nb)) ** 2 / (ca[use] + cb[use])
    ndf = int(use.sum()) - 1
    return float(x.sum()), ndf, float(C2.sf(x.sum(), max(ndf, 1)))


a, b = load(sys.argv[1]), load(sys.argv[2])
print(f"ctau(U^2=1): {a['ctau']!r} vs {b['ctau']!r} ({'identical' if a['ctau'] == b['ctau'] else 'DIFFERENT'})")
for name, k in (("channel labels", "labels"), ("charged multiplicity", "nch"), ("daughter multiplicity", "mult")):
    c, n, p = chi2_hist(a[k], b[k])
    print(f"{name}: chi2/ndf {c:.1f}/{n}, p = {p:.3f}")
for name, k in (("charged-daughter |p|", "p_ch"), ("visible momentum per decay", "pvis")):
    r = ks_2samp(a[k], b[k])
    print(f"{name}: means {a[k].mean():.4f} vs {b[k].mean():.4f} GeV, KS p = {r.pvalue:.3f}")

# decay-by-decay identity (same seed: exHad reseeds every decay, so platform differences stay local)
da, db = ({k: v for k, v in np.load(f).items()} for f in sys.argv[1:3])
na, nb = da["daughter_counts"], db["daughter_counts"]
oa, ob = np.concatenate([[0], np.cumsum(na)]), np.concatenate([[0], np.cumsum(nb)])
same = same_struct = 0
maxrel = 0.0
keys = ("px", "py", "pz", "energy")
for i in range(min(len(na), len(nb))):
    if na[i] != nb[i] or da["channel_label"][i] != db["channel_label"][i]:
        continue
    sa, sb = slice(oa[i], oa[i + 1]), slice(ob[i], ob[i + 1])
    if not np.array_equal(da["pdg"][sa], db["pdg"][sb]):
        continue
    same_struct += 1
    if all(np.array_equal(da[k][sa], db[k][sb]) for k in keys):
        same += 1
    else:
        maxrel = max(maxrel, max(float(np.max(np.abs(da[k][sa] - db[k][sb]))) for k in keys) / float(np.max(db["energy"][sb])))
print(f"decays with the same channel and daughters: {same_struct}/{len(na)}; of them bit-identical {same}, "
      f"the others differ by at most {maxrel:.1e} of the largest daughter energy")
