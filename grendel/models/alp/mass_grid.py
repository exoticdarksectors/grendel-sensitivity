"""Shared BC10 ALP mass grid with the unsupported pole points retained."""

ALP_MASS_GRID = sorted({round(x, 2) for x in (
    [0.22 + 0.02 * i for i in range(40)]
    + [1.00 + 0.05 * i for i in range(41)]
    + [1.18 + 0.01 * i for i in range(43)]
    + [3.00 + 0.10 * i for i in range(17)]
    + [3.10 + 0.02 * i for i in range(16)]
    + [4.65, 4.70, 4.75]
)})
