"""Heavy-flavour production grids: FONLL + NNPDF4.0 at pp 14 TeV.

The B-hadron and D0 double-differential cross sections the production
samplers read (``grendel/production/data/fonll/central``) are regenerated
here from FONLL v1.3.3 with the PDF set switched from the public web form's
CTEQ6.6 to ``NNPDF40_nlo_as_01180``, and validated against the public
CTEQ6.6 points. The same driver produces the scale / PDF-replica /
heavy-quark-mass variation campaign the uncertainty bands are built from.

    install    apply the patches to a FONLL source tree and build the executables
    generate   the central grids, or the variation campaign (--campaign)
    snapshot   rebuild variation_manifest.json from the grids on disk
    combine    variation envelopes (variations) and the alpha_s band (alphas)
    validate   central-reproduction gate and public CTEQ6.6 closure checks

FONLL, an LHAPDF installation and the NNPDF4.0 sets are obtained
separately (see README.md); ``GRENDEL_FONLL_DIR`` names the built tree and
``GRENDEL_FONLL_GRID_DIR`` the campaign output.
"""
