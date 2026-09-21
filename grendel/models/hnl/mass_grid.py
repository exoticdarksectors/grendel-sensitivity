"""The HNL mass grid, shared by production, channel combination and the scan.

123 points from 0.20 to 10.00 GeV, denser where the physics changes fastest:

- 0.200 - 0.500 GeV, 15 MeV steps (21 points)   kaon-threshold region
- 0.525 - 2.000 GeV, 25 MeV steps (60 points)   charm and tau thresholds
- 2.200 - 8.000 GeV, mostly 200 MeV (38 points) bottom and Bc production,
  with a 3.62 - 3.70 GeV refinement where the islands pinch shut
- 8.500 - 10.00 GeV, 500 MeV steps (4 points)   electroweak tail

Filenames encode the mass with three decimals (``grendel.io.vectors``); the
15-MeV and 25-MeV spacings alias under a two-decimal encoding.
"""
from ...io.vectors import format_mass_for_filename, parse_mass_from_filename  # noqa: F401

MASS_GRID = sorted([
    0.200, 0.215, 0.230, 0.245, 0.260, 0.275, 0.290, 0.305,
    0.320, 0.335, 0.350, 0.365, 0.380, 0.395, 0.410, 0.425,
    0.440, 0.455, 0.470, 0.485, 0.500,
    0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700,
    0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900,
    0.925, 0.950, 0.975, 1.000, 1.025, 1.050, 1.075, 1.100,
    1.125, 1.150, 1.175, 1.200, 1.225, 1.250, 1.275, 1.300,
    1.325, 1.350, 1.375, 1.400, 1.425, 1.450, 1.475, 1.500,
    1.525, 1.550, 1.575, 1.600, 1.625, 1.650, 1.675, 1.700,
    1.725, 1.750, 1.775, 1.800, 1.825, 1.850, 1.875, 1.900,
    1.925, 1.950, 1.975, 2.000,
    2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60,
    3.62, 3.63, 3.64, 3.65, 3.675, 3.685, 3.70,   # closure refinement: the islands pinch (peak_N = 3) between ~3.63 (e, mu) and ~3.69 (tau) GeV
    3.80, 4.00,
    4.20, 4.40, 4.60, 4.80, 5.00,
    5.20, 5.40, 5.50, 5.60, 5.80, 6.00,
    6.20, 6.40, 6.60, 6.80, 7.00,
    7.20, 7.40, 7.60, 7.80, 8.00,
    8.50, 9.00, 9.50, 10.00,
])

# Default per-job meson pool size; the production drivers override it.
N_EVENTS_DEFAULT = 100_000

# Full analysis range represented by the grid.
ANALYSIS_MASS_MAX = 10.0  # GeV

_labels = [format_mass_for_filename(m) for m in MASS_GRID]
assert len(set(_labels)) == len(MASS_GRID), "mass labels collide; raise the precision"
del _labels
