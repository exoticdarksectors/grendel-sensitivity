"""Sensitivity of the GRENDEL long-lived-particle detector to the PBC benchmarks.

The package is organised as one shared pipeline and one thin specification
per benchmark model:

``grendel.geometry``   the tunnel detector mesh, the four-hit reconstruction
                       and the batch ray-cast with its on-disk cache
``grendel.reco``       the acceptance Monte Carlo, the coupling scan and the
                       exclusion-band solver
``grendel.scan``       the per-mass driver every model runs through
``grendel.production`` heavy-flavour production inputs and samplers
``grendel.models``     the benchmark models (HNL, dark scalar, fermiophilic ALP)
``grendel.band``       the uncertainty-band campaigns and their combination
"""

__version__ = "0.1.0"
