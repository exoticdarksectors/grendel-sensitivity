"""Acceptance Monte Carlo, coupling scan and exclusion-band solver.

``acceptance`` builds, once per mass point, the reconstructed-decay sample
that ``scan_u2`` then reweights to every coupling; ``exclusion`` turns the
resulting yield curve into an exclusion band; ``templates`` packs
rest-frame decay templates into the bundle ``acceptance`` reads.
Importing ``acceptance`` builds the detector mesh.
"""
