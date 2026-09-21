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
| `grendel_geometry.py` | before `f2684e3` | `b16a7da27868781a714d4b0ca5cbb7a5f035314455c2bfd2e17eb0d2f1852cd2` |
| `reco_common.py` | `4a12e9774766` | `b61b38f6ccb9236873311a05022b088c76661ed0d27373e6a9d26da62eff2764` |

`grendel_geometry.py` is byte-identical to upstream. `reco_common.py` differs
in exactly one line: its sibling import `from grendel_geometry import ...`
is written `from .grendel_geometry import ...` so the two files form a
package. Nothing else is edited here; a change to the detector belongs
upstream, and refreshing means copying the file in whole, re-running every
suite and re-rendering the figures.

## Coordinate convention

`grendel_geometry` places the tunnel in the CMS frame with the beam along
`z`. Upstream later negated the beam coordinate of the tunnel vertices
(commits `f2684e3` and `6d9f9963e6ad`, `_Z_SHIFT`) to match the CMS `+z`
convention of its visualisations, describing the change as a reflection with
no effect on physics. That is true of hits, path lengths and decay
probabilities, but not of the acceptance: the tracker and scintillator
surfaces are assigned per wall (`TRACKER_SURFACES`, `SCINTILLATOR_SURFACES`),
so mirroring the tunnel without mirroring those roles sends the decay
daughters to the opposite wall and lowers the reconstructed yield by more than
an order of magnitude. This package therefore keeps the coordinates every
result was computed with, and does not adopt the negation. Anything that maps
decays into another frame (a Geant4 model of the gallery, say) has to fix the
parity between the two descriptions explicitly.
