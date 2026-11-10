#!/usr/bin/env python3
"""
hnl/config_mass_grid.py

HNL mass grid used by the production drivers.

The grid (116 points, 0.20 - 10.00 GeV) is denser where the physics changes
fastest:

- 0.200 - 0.500 GeV step 15 MeV    21 points   kaon-threshold region
- 0.525 - 2.000 GeV step 25 MeV    60 points   charm and tau thresholds
- 2.200 - 8.000 GeV mostly 200 MeV 31 points   bottom and Bc production
- 8.500 - 10.00 GeV step 500 MeV    4 points   electroweak tail

The same grid is used by production, channel combination, and analysis.

`format_mass_for_filename(m)` controls the filename encoding:
`mN_{m:.3f}.csv` with the decimal replaced by `p` (e.g. 1.025 -> mN_1p025).
Three decimals are required: the grid has 15-MeV and 25-MeV spacings, so a
two-decimal encoding aliased neighbouring points and drifted the label up to
5 MeV from the true `m_N`. `parse_mass_from_filename` is the exact inverse.
"""

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
    2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
    4.20, 4.40, 4.60, 4.80, 5.00,
    5.20, 5.40, 5.50, 5.60, 5.80, 6.00,
    6.20, 6.40, 6.60, 6.80, 7.00,
    7.20, 7.40, 7.60, 7.80, 8.00,
    8.50, 9.00, 9.50, 10.00,
])

# Default per-job meson pool size. Overridden by run_all.py / driver CLI flags.
N_EVENTS_DEFAULT = 100_000

# Hard cap on emitted HNL events per (flavor, channel, mass). 0 = unlimited.
# Reserved for future filtering; currently unused by the drivers.
MAX_SIGNAL_EVENTS = 0


def format_mass_for_filename(mass):
    """Filename encoding: 1.025 -> '1p025' (three decimals, dot as 'p')."""
    return f"{mass:.3f}".replace('.', 'p')


def parse_mass_from_filename(label):
    """Inverse of ``format_mass_for_filename``.

    Accepts either a bare label (``'1p025'``) or a full stem
    (``'mN_1p025'``) and returns the mass in GeV (``1.025``).
    """
    stem = label[3:] if label.startswith("mN_") else label
    return float(stem.replace('p', '.'))


# Fail fast if the grid ever grows a point whose label collides with another:
# a collision would silently overwrite a CSV. Three decimals (1-MeV resolution)
# safely separate the 15-MeV/25-MeV grid spacings used here.
_labels = [format_mass_for_filename(m) for m in MASS_GRID]
assert len(set(_labels)) == len(MASS_GRID), (
    "format_mass_for_filename produced colliding labels for MASS_GRID; "
    "increase the decimal precision before adding the new point(s)."
)
del _labels


if __name__ == "__main__":
    print(f"MASS_GRID: {len(MASS_GRID)} points ({min(MASS_GRID):.2f} - {max(MASS_GRID):.2f} GeV)")
    print(f"  Below 2 GeV: {len([m for m in MASS_GRID if m < 2])} points")
    print(f"  Above 2 GeV: {len([m for m in MASS_GRID if m >= 2])} points")
