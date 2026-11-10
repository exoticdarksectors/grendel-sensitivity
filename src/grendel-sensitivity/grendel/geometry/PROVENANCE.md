# Detector geometry and reconstruction

`grendel_geometry.py` and `reco_common.py` are the GRENDEL detector's own
description -- the tunnel mesh, the fiducial volume, the surface
classification and the bounded four-hit vertex reconstruction with its
selection -- shared with the detector-optimisation and cosmic-background
studies in [exoticdarksectors/llpatcolliders](https://github.com/exoticdarksectors/llpatcolliders)
(`higgs/`). They are taken from there whole, so signal and background are
reconstructed by one definition.

| File | Upstream commit | sha256 of the upstream file |
|---|---|---|
| `grendel_geometry.py` | `7b12c7e59a1d` | `b16a7da27868781a714d4b0ca5cbb7a5f035314455c2bfd2e17eb0d2f1852cd2` |
| `reco_common.py` | `4a12e977476b` | `b61b38f6ccb9236873311a05022b088c76661ed0d27373e6a9d26da62eff2764` |

`grendel_geometry.py` was byte-identical to upstream up to release v2.2. Since
October 2026 its tunnel is UP542 as built, read from the mirror in `data/up542/`
(`../../GEOMETRY_UP542.md`): measured section, as-built centreline and floor,
official CMS frame. The old centreline was a mirror image of the gallery; the
as-built one is walked from the south end towards PX56, so that the wall roles,
defined left and right of the walking direction, stay those of the published
design (`TRACKER_SURFACES` = ceiling + left wall, `SCINTILLATOR_SURFACES` = floor +
right wall, the right wall facing IP5 on 80 % of the gallery). Walked from PX56,
the trackers would face IP5 and 12 % of the daughters leaving the volume along
the LLP direction would end on the veto instead. `reco_common.py` differs
only in its imports: the sibling import `from grendel_geometry import ...` is
written `from .grendel_geometry import ...` so the two files form a package,
and `points_in_fiducial` comes from `fiducial.py` instead of
`grendel_geometry`. Upstream's `points_in_fiducial` tests the cross-section
of the nearest centreline segment, while the mesh places each ring
perpendicular to the node tangent; at the bends (up to 31 deg) it accepts
0.4 % of volume outside the mesh, up to 20 cm deep and never the reverse,
which biased the reconstructed yield up by 0.2-0.7 %. `fiducial.py` tests
containment in the mesh itself (vertical-ray parity against its triangles,
identical to `mesh.contains`), so the ray-cast and the vertex and hit
selections share one volume. Nothing else is edited in `reco_common.py`; a change to the
detector belongs upstream (this one should go there too), and refreshing
means copying the file in whole, re-running every suite and re-rendering
the figures.

## Coordinate convention

Since October 2026 the tunnel is in the official CMS physicist frame with the
correct handedness (`GEOMETRY_UP542.md`). Up to v2.2 it was the mirror image of
the gallery, with the beam coordinate opposite in sign to CMS `+z`. An upstream
branch (commits `f2684e3` and `6d9f9963e6ad`, `_Z_SHIFT`; not on upstream
`main`) negated the beam coordinate of the old vertices, describing the change
as a reflection with no effect on physics. That is true of hits, path lengths
and decay probabilities, but not of the acceptance: the tracker and
scintillator surfaces are assigned per wall (`TRACKER_SURFACES`,
`SCINTILLATOR_SURFACES`), so mirroring the tunnel without mirroring those roles
sends the daughters that reach a side wall (about a fifth) to the opposite one
(on the BC4, BC10 and HNL inputs: peak yield 1-9 % lower, lower coupling edge
2-8 % higher). The as-built tunnel keeps the roles by its walking direction
(above). Anything that maps decays into another frame (a Geant4 model of the
gallery, say) has to fix the parity between the two descriptions explicitly.
