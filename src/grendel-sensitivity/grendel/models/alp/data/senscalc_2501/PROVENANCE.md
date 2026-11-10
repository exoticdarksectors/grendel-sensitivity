# SensCalc arXiv:2501.04525 decay inputs

This directory is the reviewable export of the universally fermion-coupled
ALP decay data used by GRENDEL BC10. The source is the official SensCalc
release `v1.4`, pinned at commit
`4f7210b6e566f5510e2c7269adec3744a6369880` (2026-10-11).

Upstream repository: <https://github.com/maksymovchynnikov/SensCalc>

SensCalc methodology: arXiv:2305.13383.

Decay phenomenology: arXiv:2501.04525. This description is invariant under
the chiral-rotation parameter choice and includes mixing with heavy
pseudoscalar excitations. SensCalc v1.4 ships it as the `Default` decay
description (the files named `2501.04525` in v.1.3.3).

Against the v.1.3.3 export used until 2026-10-13: the mu mu width is
identical, the total width differs by at most 0.41 % (at 2.17 GeV), and the
a -> K_L K_S pi0 channel (branching ratio at most 0.4 %, 1.13-2.17 GeV) is
zero in v1.4. The squared matrix elements are unchanged.

## Source files

The binary Wolfram MX dumps are not copied into this repository. The exporter
checks these exact SHA-256 values before decoding them:

| Input | SHA-256 |
|---|---|
| `Widths-model-ALP-fermion-scale-1000.-GeV-Default.m` | `8e438677df92ed3354e7004a667604bd665b182d0ba2780f0ec2262c5caf96d1` |
| `Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-Default.m` | `def3041c57db36f6074deab0d373a4fb893453649ac0fca43ec45b2fc6d4e484` |
| `Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-Default.m` | `f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed` |
| `codes/Acceptances/ALP-fermion.nb` | `d44282084ffe2565c5811254e33265f2575ac80f4bc4b168ff841fda491d8a66` |

All three inputs are retained together. Widths determine the lifetime and
branching ratios determine the visible-channel mixture. The squared matrix
elements supply the primary three-body shapes. The publication generator in
`templates_pythia.py` uses the exact exclusive branching
mixture, reweights flat three-body phase space with those squared matrix
elements, and lets Pythia decay unstable daughters. The light-quark surrogate
for the two-gluon mode remains an explicitly documented approximation.

The exporter retains every upstream display label and expression, but it also
assigns canonical machine names to the unambiguous leptonic, diphoton, and
total-width entries. Model code must use those canonical fields rather than
parse Mathematica box notation or localized Greek characters.

The pinned dumps record the Wolfram system IDs `MacOSX-ARM64` (widths,
branching ratios) and `Windows-x86-64` (matrix elements) in their MX headers. MX is documented as system-dependent, so the exporter reports the
source system ID and records the decoder `$SystemID` in the
`EXPORT_MANIFEST.json` it writes (not committed here, nor is
`widths_raw_senscalc.csv`). The committed bundle was successfully decoded and
validated with Wolfram Engine 15.0 on `MacOSX-ARM64`; this demonstrated import
is recorded in the manifest rather than inferred from the file header. See the
[Wolfram MX format documentation](https://reference.wolfram.com/language/ref/format/MX.html).

## Regeneration

Obtain and activate the free Wolfram Engine, then run:

```bash
git clone --branch v1.4 --depth 1 \
  https://github.com/maksymovchynnikov/SensCalc /tmp/SensCalc-v1.4
wolframscript -activate
python -m grendel.models.alp.tools.export_senscalc_2501 /tmp/SensCalc-v1.4
```

On the Homebrew macOS installation, the wrapper automatically selects the
nested kernel at:

```text
/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel
```

Use `--check-only` to validate the SensCalc commit and source hashes without a
Wolfram license. It also prints the system IDs embedded in the MX headers.

## Coupling conversion

SensCalc applies its raw width coefficient as

```text
Gamma = (g_Y / (2 v_h))^2 * coefficient_raw.
```

The GRENDEL BNT-axis convention is `1/f_BNT = g_Y/v_h`, so the exported BNT
coefficient is `coefficient_raw/4`. The branching ratios are convention
independent. The tables are the default decay-model input. Validation checks
the raw-to-BNT factor row by row and the exact `BR(a -> mu mu) = 0.2022433833`
anchor at 1 GeV. The older approximately 9% anchor came from the superseded
GKOZ/ALPINIST digitization, read without ALPINIST's 1/4, and must not be
imposed on the 2501 tables.

SensCalc's charged/no-ECAL channel policy is reproduced from
`codes/EventCalc/DecayProductsSampler.nb`: the diphoton, `3 pi0`, and
`2 K_L pi0` channels are excluded from the visible branching fraction. The
remaining exclusive branching ratios sum to the visible fraction that only the
legacy proxy templates use; the publication templates keep the full branching
and leave visibility to the GRENDEL reconstruction.

The upstream MX list repeats the charged and neutral `K* K*` processes at
positions 11/29 and 12/24. Each pair has identical products and a byte-identical
branching expression. SensCalc installs these entries as Mathematica
`DownValues` keyed by process name, so the repeated definition overwrites the
first rather than contributing a second branching fraction. The positional CSV
retains all source rows for auditability; `model.py` excludes channels 24 and
29 from sums to reproduce the SensCalc semantics. A full-table test enforces
that the resulting visible fraction stays in `[0, 1]`.

SensCalc is used as a verified
physics-input source; its detector acceptance is not substituted for the
GRENDEL geometry and reconstruction chain.
