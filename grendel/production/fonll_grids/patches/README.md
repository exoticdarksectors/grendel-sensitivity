# Patches to FONLL 1.3.3

Unified diffs against FONLL 1.3.3 (`alisw/fonll` at `c7086e49`); `PINS.json`
holds the SHA-256 of each file before and after patching.
`python -m grendel.production.fonll_grids.install` applies them (by hand:
`patch -p1 < <diff>` from the FONLL root). They widen array capacities and
output formats, add the BCFY fragmentation functions (options 4, 5 and 8;
hep-ph/9409316) and adapt the build to gfortran and LHAPDF.
