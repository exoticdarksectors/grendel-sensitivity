"""Shared BC10 ALP mass grid with the unsupported pole points retained.

The sensitivity output needs marker rows at the eta and eta-prime poles so
plotting cannot bridge them. Generators inspect ``model`` (or the standalone
window constants in the ROOT-only template generator) and skip those points.
"""

ALP_MASS_GRID = sorted({round(x, 2) for x in (
    [0.22 + 0.02 * i for i in range(40)]
    + [1.00 + 0.05 * i for i in range(41)]
    # Resolve the heavy-pseudoscalar width structure and narrow sensitivity
    # windows without interpolating between independently simulated masses.
    + [1.18 + 0.01 * i for i in range(43)]
    + [3.00 + 0.10 * i for i in range(17)]
    # Resolve the N_signal=3 closure of the high-mass island.
    + [3.10 + 0.02 * i for i in range(16)]
    + [4.65, 4.70, 4.75]
)})
