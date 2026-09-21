# ANUBIS — original-source provenance

**Latest HNL projection checked:** arXiv:2606.26862v1, posted 25 Jun 2026,
"Projected sensitivity of the ANUBIS detector to heavy neutral leptons".

**Older HNL projections:** arXiv:2001.04750 and arXiv:2512.14942.

**Hadronic/meson-production studies:**

- Hirsch and Wang, arXiv:2001.04750, DOI: 10.1103/PhysRevD.101.055034.
- Wang and Zhang, arXiv:2512.13011, DOI: 10.1103/fyd4-gk83.

**Figure containing the latest contours:** Fig. 4 of arXiv:2606.26862v1.
The source archive is unpacked at `sources/arXiv-2606.26862v1/`; the standalone
Matplotlib figure PDFs used for extraction are:

- `sources/arXiv-2606.26862v1/img/HNL_electron_ANUBIS_Combination_total_sigEff0p5.pdf`
- `sources/arXiv-2606.26862v1/img/HNL_muon_ANUBIS_Combination_total_sigEff0p5.pdf`

**Extracted data:**

- `ANUBIS_2026_Ue.dat`: BC6/electron-coupled Majorana HNL, solid red total contour.
- `ANUBIS_2026_Umu.dat`: BC7/muon-coupled Majorana HNL, solid red total contour.
- `ANUBIS_WangZhang_2025_ceiling_Ue.dat`: ANUBIS-ceiling hadronic/meson
  production contour from arXiv:2512.13011, Fig. 1, solid dark-blue curve.

**Important scope note:** arXiv:2606.26862v1 studies electron and muon single-flavor
couplings only. The tau-only BC8 scenario is explicitly reserved for future study,
so no latest ANUBIS tau contour is included here.

**Hadronic/meson scope note:** arXiv:2512.13011 studies HNLs mixed with the
electron neutrino only and includes charm/bottom meson production only. It notes
that muon-only sensitivity should be similar apart from threshold effects, but no
muon curve is plotted. For the muon plot in this repository, the available
`curves/data/hnl/vector/ANUBIS_Umu.dat` curve is a PBC BC7 vector
digitization, not a direct Wang-Zhang source curve.

**Method:** vector extraction from the standalone Matplotlib PDF path operators,
with log-axis calibration from the figure grid. No public machine-readable curve
data was found in the arXiv source archive.

The Wang-Zhang source archive is unpacked at `sources/arXiv-2512.13011v1/`; the
figure used for the electron hadronic/meson curve is
`sources/arXiv-2512.13011v1/8_sensitivity_VvSq_vs_mass_DandB_merged_LHCandSHiP.pdf`.

**Cached PDFs:**

- 2001.04750: `<work dir>/curves/original_sources_pdfs/2001.04750.pdf`
- 2512.14942: `<work dir>/curves/original_sources_pdfs/2512.14942.pdf`
- 2606.26862: `<work dir>/curves/original_sources_pdfs/2606.26862.pdf`
