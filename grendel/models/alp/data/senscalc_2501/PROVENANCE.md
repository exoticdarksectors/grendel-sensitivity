# SensCalc arXiv:2501.04525 decay inputs

This directory is the reviewable export of the universally fermion-coupled
ALP decay data used by GRENDEL BC10. The source is the official SensCalc
release `v.1.3.3`, pinned at commit
`0bca050633aae16e148d47f21840fa07ff4b8724`.

Upstream repository: <https://github.com/maksymovchynnikov/SensCalc>

SensCalc methodology: arXiv:2305.13383.

Decay phenomenology: arXiv:2501.04525. This description is invariant under
the chiral-rotation parameter choice and includes mixing with heavy
pseudoscalar excitations. SensCalc `v.1.3.3` is used because its release notes
specifically record a fix to the decay-mode description of universally
fermion-coupled ALPs.

## Source files

The binary Wolfram MX dumps are not copied into this repository. The exporter
checks these exact SHA-256 values before decoding them:

| Input | SHA-256 |
|---|---|
| `Widths-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `d12fb78d28edff0aa081c9fb66d829b42b4ec71202684019d7b9047ecb40b869` |
| `Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `34a09eed87d081bfffe79b464094741454022f79478c5b28bc236bc361049286` |
| `Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed` |
| `codes/Acceptances/ALP-fermion.nb` | `8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4` |

All three inputs are retained together. Widths determine the lifetime and
branching ratios determine the visible-channel mixture. The squared matrix
elements supply the primary three-body shapes. The publication generator in
`generate_decay_templates_pythia.py` uses the exact exclusive branching
mixture, reweights flat three-body phase space with those squared matrix
elements, and lets Pythia decay unstable daughters. The light-quark surrogate
for the two-gluon mode remains an explicitly documented approximation.

The exporter retains every upstream display label and expression, but it also
assigns canonical machine names to the unambiguous leptonic, diphoton, and
total-width entries. Model code must use those canonical fields rather than
parse Mathematica box notation or localized Greek characters.

The pinned dumps record the Wolfram system ID `Windows-x86-64` in their MX
headers. MX is documented as system-dependent, so the exporter reports the
source system ID and records the decoder `$SystemID` in
`EXPORT_MANIFEST.json`. The committed bundle was successfully decoded and
validated with Wolfram Engine 15.0 on `MacOSX-ARM64`; this demonstrated import
is recorded in the manifest rather than inferred from the file header. See the
[Wolfram MX format documentation](https://reference.wolfram.com/language/ref/format/MX.html).

## Regeneration

Obtain and activate the free Wolfram Engine, then run:

```bash
git clone --branch v.1.3.3 --depth 1 \
  https://github.com/maksymovchynnikov/SensCalc /tmp/SensCalc-v1.3.3
wolframscript -activate
python alp_fermion/tools/export_senscalc_2501.py /tmp/SensCalc-v1.3.3
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
the raw-to-BNT factor row by row, the exact `BR(a -> mu mu) = 0.2022433833`
anchor at 1 GeV, and the charm-region structure near 2.58 GeV. The older
approximately 9% anchor came from the superseded GKOZ/ALPINIST digitization and
must not be imposed on the 2501 tables.

SensCalc's charged/no-ECAL channel policy is reproduced from
`codes/EventCalc/DecayProductsSampler.nb`: the diphoton, `3 pi0`, and
`2 K_L pi0` channels are excluded from the visible branching fraction. The
remaining exclusive branching ratios are summed before the separate GRENDEL
decay-template and detector-acceptance calculation.

The upstream MX list repeats the charged and neutral `K* K*` processes at
positions 11/29 and 12/24. Each pair has identical products and a byte-identical
branching expression. SensCalc installs these entries as Mathematica
`DownValues` keyed by process name, so the repeated definition overwrites the
first rather than contributing a second branching fraction. The positional CSV
retains all source rows for auditability; `model.py` excludes channels 24 and
29 from sums to reproduce the SensCalc semantics. A full-table test enforces
that the resulting visible fraction stays in `[0, 1]`.

The production-channel comparison is tracked separately in
[`PRODUCTION_AUDIT.md`](PRODUCTION_AUDIT.md). SensCalc is used as a verified
physics-input source; its detector acceptance is not substituted for the
GRENDEL geometry and reconstruction chain.
