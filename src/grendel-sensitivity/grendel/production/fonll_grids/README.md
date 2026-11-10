# FONLL meson grids

Regeneration of the B-hadron and D0 double-differential production cross
sections at pp 14 TeV with FONLL v1.3.3 and `NNPDF40_nlo_as_01180`, and
the scale / PDF-replica / heavy-quark-mass variation campaign behind the
uncertainty bands. The public FONLL web form offers no NNPDF4.0 set; this
package reproduces its conventions with NNPDF4.0 and validates against the
public CTEQ6.6 points. The two central grids it produced are tracked in
`../data/fonll/central` (provenance in `../data/fonll/PROVENANCE.md`); the
campaign grids are regenerable and large, and are not.

## Central setup

| Quantity | Value |
| --- | --- |
| Process | pp, sqrt(s) = 14 TeV, ebeam1 = ebeam2 = 7000 GeV |
| Calculation | FONLL v1.3.3, NLO + NLL |
| PDF | `NNPDF40_nlo_as_01180`, LHAPDF id 331700 |
| Heavy-quark masses | m_b = 4.75 GeV, m_c = 1.50 GeV |
| Scales | mu_R = mu_F = sqrt(m^2 + pT^2) |
| Grid | 100 pT in [0, 50] GeV x 100 y in [-3, 3]; `pT y dsigma/dpT/dy` in pb/GeV |
| Bottom fragmentation | Kartvelishvili `(1-z) z^alpha`, alpha = 24.2 (the public B-hadron default with `n5moment`) |
| Charm fragmentation | BCFY pseudoscalar, r = 0.1, plus a calibrated D* -> D0 feeddown; the public D0 convention |
| Fragmentation frame | `ifrframe = 1` (y = 0 frame, the public default for `dsigma/dpT/dy`) |

## Setup

1. Python >= 3.11 with numpy (the generator itself is standard library).
2. `gfortran`, GNU make, and an LHAPDF >= 6 whose `lhapdf-config` is on
   `PATH`.
3. The FONLL v1.3.3 source tree: the `alisw/fonll` mirror at the commit in
   `patches/PINS.json` (`python third_party/fetch.py fonll`; the authors'
   page is <https://www.lpthe.jussieu.fr/~cacciari/fonll/>). Put it
   at `third_party/fonll` or point `GRENDEL_FONLL_DIR` at it.
4. The PDF sets: `lhapdf install NNPDF40_nlo_as_01180`; `cteq66` for the
   charm feeddown calibration and the validation runs;
   `NNPDF40_nlo_as_01170` and `_01190` for the alpha_s band.
   `GRENDEL_LHAPDF_DATA`, when set, is exported as `LHAPDF_DATA_PATH`.

Then apply the patches and build the two executables:

```
python -m grendel.production.fonll_grids.install --build
```

On Linux the installer passes the Darwin compiler line,
`F77=gfortran -Wall -fno-automatic -fallow-argument-mismatch`, and
`SYSOBJ=` to make, so the build matches the macOS one the grids were
produced with: the patched Makefile moves only its Darwin block to
gfortran, and its Linux block would call `g77` and link the x86-only
floating-point-trap object `trapfpe.o`.

`patches/` holds the four unified diffs against the pristine v1.3.3 files
with their before/after SHA-256; `install` refuses a tree that is neither
pristine nor already patched. `patches/README.md` says what each hunk does.

## Central grids

```
python -m grendel.production.fonll_grids.generate --pdf nlo --quark bottom --grid-workers 4
python -m grendel.production.fonll_grids.generate --pdf nlo --quark charm  --grid-workers 4
```

Grids land in `--out-dir` (`GRENDEL_FONLL_GRID_DIR`, default
`<work>/fonll_grids/output`) with run scratch and logs beside it.
`--grid-workers N` splits the rapidity axis into `N` contiguous chunks
computed by independent FONLL processes and re-stitched before
fragmentation; there is no shared state across chunks.

Charm is two `fragmfonll` passes over one quark grid, the direct BCFY
pseudoscalar stream (mode 5) and the BCFY vector stream (mode 4) evaluated
at the D* parent kinematics, combined with a single global weight (the
mix is renormalised to fragmentation fraction 1) fitted to nine public
FONLL CTEQ6.6 D0 points at y = 0 and cached in
`charm_feeddown_calibration.json`. The weight is a property of the
fragmentation model, so every charm variation reuses the one calibration.
`generate --calibrate-only` fits and caches it alone; run it once before
starting charm grids in parallel, which would otherwise all fit it in the
same run directory.

Before trusting a fresh build, confirm it reproduces the tracked centrals:

```
python -m grendel.production.fonll_grids.validate gate --quark bottom --grid-workers 6
python -m grendel.production.fonll_grids.validate public-points   # needs generate --pdf cteq66 --quark bottom
```

## Variation campaign

```
python -m grendel.production.fonll_grids.generate --campaign --max-parallel 6 --grid-workers 1 --compress-logs
python -m grendel.production.fonll_grids.combine variations
```

The bare `--campaign` runs the full default set per quark: 7-point scale
(`(muR, muF)` over {0.5, 1, 2} without the antipodal extremes), NNPDF4.0
replica members 1-100 (LHAPDF id `331700 + member`), and
`m_b = 4.75 +/- 0.25`, `m_c = 1.5 -/+ 0.2` GeV; 218 grids in all. Narrow
it with `--scale-variations`, `--mass-variations` and `--pdf-members 0-30`.
Each grid is about two single-core hours, so size `--max-parallel` to the
core count and keep `--grid-workers 1`. The run is resumable with
`--reuse-existing-grids` (a raw quark grid an interrupted run left
incomplete is rebuilt). FONLL's grid driver reads the two scale factors
as `(ffact = muF, fren = muR)`; the generator maps each point accordingly.

Every grid records its variation in the filename and header (`ffact`,
`fren`, `lhapdf_id`, `lhapdf_member`, `heavy_quark_mass_GeV`,
`variation_kind`, `variation_tag`). `variation_manifest.json` records the
FONLL revision, patches, PDF set and members, scale points, masses, grid
bounds, the charm calibration, and a SHA-256 and integrated cross section
per grid. The model band campaigns (`grendel.band.*`) read that manifest
and verify the hashes of the grids they use. `snapshot` rebuilds the
manifest from the grids on disk so a partial campaign can be combined and
banked at any time.

`combine variations` writes envelope grids under `envelopes/`: `scaleup`/
`scaledn` (max/min over the scale set), `pdfmean`, `pdfup`/`pdfdn`
(central +/- replica std, ddof = 1, member 0 excluded), `massup`/`massdn`,
and `combup`/`combdn` (the three in uncorrelated quadrature). Bands are
only emitted for axes with at least two grids.

The alpha_s band uses the `NNPDF40_nlo_as_01170` / `_01190` companions
through the plain path, then the PDF4LHC half-difference:

```
python -m grendel.production.fonll_grids.generate --pdf nlo_as_01170 --pdf nlo_as_01190 --quark bottom --quark charm
python -m grendel.production.fonll_grids.combine alphas
```

`campaign.sh` chains the whole thing for an unattended run: its Phase 0
fits the charm calibration once, serially (`--calibrate-only`), before the
parallel phases, and every replica stage banks a complete manifest and
envelopes.

**Using the bands.** The `combined` envelope is scale + PDF + mass and
does not include alpha_s; add the `alphas` band in quadrature for the
total. The grids are per-quark rates, `(q + qbar)/2`, at fragmentation
fraction 1 (bottom: any B hadron; charm: the D0-convention shape) — a
species rate is `2 * sigma * f * BR` with the measured fraction. Integrated
charm sigma *decreases* with alpha_s, pointwise only below pT ~ 2.5 GeV and
above ~ 22 GeV (the alpha_s–gluon anti-correlation of the NNPDF4.0 fit at
the ~1.5 GeV charm scale outweighs the explicit alpha_s growth). Envelopes
carry an `envelope_band` header and are pointwise diagnostics: the
production sampler refuses them, and only coherent individual grids may be
sampled.

## Citing

Cite FONLL as the public web form asks — Cacciari, Greco, Nason, JHEP 9805
(1998) 007 [hep-ph/9803400]; Cacciari, Frixione, Nason, JHEP 0103 (2001)
006 [hep-ph/0102134] — and the NNPDF4.0 release, Eur. Phys. J. C 82 (2022)
428 [arXiv:2109.02653]. The BCFY fragmentation functions are Braaten,
Cheung, Fleming, Yuan, Phys. Rev. D51 (1995) 4819 [hep-ph/9409316].
