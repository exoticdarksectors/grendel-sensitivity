# Provenance of the central grids

    central/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat
    central/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat

## How they were generated

- Code: FONLL v1.3.3 (April 2014), the `alisw/fonll` mirror at commit
  `c7086e49141cf6705cf7a4bc5f7d0b3a38673203`, with the modifications in
  `grendel/production/fonll_grids/patches`.
- PDF backend: LHAPDF 6.5.6.
- Driver: `grendel.production.fonll_grids.generate --pdf nlo --quark {bottom,charm}`.

## Bottom grid

Public FONLL B-hadron default fragmentation, Kartvelishvili `(1-z) z^alpha`
with central `alpha = 24.2`, fragmentation fraction 1. A single direct
production stream from the patched `fragmfonll` on the 100 x 100 (pT, y)
grid, internal quark-side fragmentation cutoffs at FONLL defaults.

Closure against the public FONLL v1.3.2 CTEQ6.6 web form, queried
2026-05-29:

    B hadron, pT = 5 GeV, y = 0
        local CTEQ6.6 run :  7.535277e6 pb/GeV
        public CTEQ6.6    :  7.532900e6 pb/GeV
        relative          : +0.032%

The local CTEQ6.6 run uses the same code path as the NNPDF4.0 grid; only
the LHAPDF set differs, so the 3e-4 agreement at the reference point bounds
the systematic introduced by the patched code itself.

## Charm grid

The public FONLL D0 convention, the sum of two contributions computed
separately by the patched code:

1. direct BCFY pseudoscalar (mode 5), r = 0.1, fragmentation fraction 1;
2. direct BCFY vector (mode 4), r = 0.1, fragmentation fraction 1,
   convolved with a collinear two-body `D* -> D0` decay using PDG masses and
   BR(D*0 -> D0 pi0), BR(D*0 -> D0 gamma), BR(D*+ -> D0 pi+).

The two streams are combined with one global weight,
`charm_feeddown_vector_to_direct_weight = 1.21525336412`, fitted to nine
public FONLL CTEQ6.6 D0 points at y = 0 and pT = 1, 5, 8, 15, 22, 29, 36,
43, 50 GeV; the maximum absolute relative residual over them is 1.94e-3.

Closure on the direct streams alone, for reference:

    D* vector,  pT = 5 GeV, y = 0, r = 0.1
        local CTEQ6.6 run :  4.0129e7 pb/GeV
        public CTEQ6.6    :  4.0129e7 pb/GeV
        relative          :  match at the displayed precision

    D0 direct pseudoscalar, pT = 5 GeV, y = 0, r = 0.1 (no feeddown)
        local CTEQ6.6 run :  3.679271e7 pb/GeV
        public CTEQ6.6    :  3.522100e7 pb/GeV
        relative          : +4.46%

The +4.46 % offset of the direct stream is what the calibrated feeddown
absorbs. Without the feeddown convolution the grid is not consistent with
the public D0 convention and should not be used.

## Conventions

- The grid is `pT y dsigma/dpT/dy` in pb/GeV. The public web form also
  offers `dsigma/dpT^2/dy` in pb/GeV^2, which differs by `1/(2 pT)`.
- FONLL returns `(q + qbar)/2`. Studies that count produced mesons multiply
  by 2.
- Fragmentation frame `ifrframe = 1` (y = 0 frame). The meson mass `xmh`
  enters only the `szmin` kinematic guard in `fragmfonll`,
  `(pT^2 + xmh^2) cosh(y)^2 > sh/4`. Across this grid (`pT <= 50 GeV`,
  `|y| <= 3`, `sh/4 = 4.9e7 GeV^2`) the left-hand side is at most about
  `2.5e5 GeV^2`, two orders of magnitude below the threshold, so the grid
  is insensitive to `xmh` here: regenerating with `xmh = M_B0 / M_D0` gives
  byte-identical numerics, and the samplers put the meson on shell with
  its PDG mass at sampling time. Extending the grid toward the kinematic
  edge would re-activate the guard and require re-verification.
- The `internal_quark_grid` header names the raw quark grid the meson grid
  was fragmented from; `fragmentation_fraction: 1` is a convention label,
  not a literal f = 1 (see below).

## Charm species

The charm grid is a *D0-convention* shape calibrated to the public FONLL
D0 output. Total charm-meson rates downstream still use measured species
fractions f(D0), f(D+), f(Ds), f(Lambda_c); the grid provides the pT–y
shape under the D0 convention only.

## Uncertainty bands

The theory band is produced by the campaign tooling in
`grendel/production/fonll_grids` and is not tracked (the grids are
regenerable and large). Each product records a SHA-256 per grid in its
manifest: `generate --campaign` (218 coherent grids per full run:
central + 7-point scale + 100 NNPDF4.0 replicas + 2 heavy-quark-mass
points per quark, catalogued in `variation_manifest.json`),
`combine variations` (scale / PDF / mass / combined envelopes under
`envelopes/`), `combine alphas` (the strong-coupling band from the
`NNPDF40_nlo_as_01170` / `_01190` companion centrals).

Caveats:

1. **alpha_s is a separate band.** The `combined` envelope is scale + PDF
   + mass in quadrature only; add the `alphas` band in quadrature for the
   total.
2. **The grids are absolute meson cross sections.** Bottom is the
   inclusive B-hadron rate; charm is the D0-convention rate with the
   c -> D0 normalisation absorbed into the calibrated feeddown weight. Do
   not re-apply a fragmentation fraction f(b -> B) or f(c -> D0). The
   species fractions above apply only when *summing* D0/D+/Ds/Lambda_c.
3. **Charm sigma decreases with alpha_s** (sigma(0.117) > sigma(0.118) >
   sigma(0.119)), opposite to bottom. Physical: at the ~1.5 GeV charm scale
   the anti-correlation between alpha_s and the gluon in the NNPDF4.0 fit
   outweighs the explicit alpha_s growth (verified against the LHAPDF
   `.info` metadata, no set mislabelling). Integrated band: +/-1.49 %
   (bottom), +/-1.88 % (charm) of the central.

## The patches

Documented hunk by hunk in `grendel/production/fonll_grids/patches/README.md`:
the rapidity capacity of the grid driver and the fragmentation code
raised to 120 nodes so a 100-node rapidity grid plus its sentinel fits;
the BCFY vector and pseudoscalar fragmentation functions (modes 4 and
5/8) after Braaten, Cheung, Fleming, Yuan, Phys. Rev. D51 (1995) 4819;
higher-precision output formats; and the LHAPDF link target in the
Makefiles. No physics formula of FONLL is altered.
