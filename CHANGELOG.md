# Changelog

## 1.0 — 2026-09-21

First release: the sensitivity machinery of arXiv:2609.00152.

- Detector geometry, two-track vertex reconstruction and the acceptance
  Monte Carlo shared by every model.
- One scan driver (`grendel.scan`) with per-model numerical policies for the
  heavy neutral lepton (BC6–BC8), the dark scalar (BC4) and the fermiophilic
  ALP (BC10); islands solved at any number of signal thresholds.
- Production: FONLL heavy-flavour grids and their generator, HNL channels
  (mesons, baryons, kaons, tau, W/Z), B → K S and B → K a.
- Decay templates from exHad, FairShip and the analytic and proxy backends.
- Uncertainty campaigns (FONLL scale, PDF replicas and heavy-quark masses;
  decay-model, B_c and numerical-control nuisances) and their combination.
- The comparison-curve renderer with the digitised reference curves and
  the tools that produced them.
