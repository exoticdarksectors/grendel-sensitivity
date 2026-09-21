# Patches to FONLL v1.3.3

Unified diffs against the pristine FONLL 1.3.3 import (the `alisw/fonll`
mirror at commit `c7086e49`, identical to the authors' tarball). The
SHA-256 of each file before and after patching is in `PINS.json`;
`python -m grendel.production.fonll_grids.install` applies them and
checks both. Apply by hand from the FONLL root with `patch -p1 < <diff>`.

| Diff | File | What changes |
|---|---|---|
| `misc1_fonllgrid.f.diff` | grid driver | rapidity array capacity `nymax` 100 → 120; the `.out` and `.outlog` files are opened in append mode instead of being wound to the end with `toend`; the raw grid record is written with `es18.10` instead of `d12.6` |
| `misc1_fragmfonll.f.diff` | fragmentation convolution | rapidity capacity `nymx` 50 → 120 at the four declarations that share the COMMON layout; BCFY fragmentation options 4 (vector) and 5/8 (pseudoscalar); an explicit stop on an unknown option; the output records written with `es18.10` instead of `d10.4`; the free-format number reader widened from 15 to 24 characters to read them back |
| `misc1_Makefile.diff` | build | `gfortran` on Darwin (`-fallow-argument-mismatch`); the `fonllgridlha` target linking the grid driver against LHAPDF through `lhapdf-config`; the objects it needs (`lhapdfif.o`, `eks98.o`, `eps09.o`, `struv_plus_avgo.o`) |
| `main_Makefile.diff` | build | `-Wl,-rpath` to the LHAPDF library directory |

No physics formula is altered by the capacity or format changes: the
arrays are storage, the loop bounds are driven by the number of nodes read
from the input card, and the wider output format only stops the raw grid
being truncated to six digits on its way into the fragmentation step.

## The BCFY fragmentation functions

`fragmfonll.f::fragfun` gains the S-wave heavy–light fragmentation
functions of Braaten, Cheung, Fleming and Yuan, Phys. Rev. D51 (1995)
4819 [hep-ph/9408231], with the parameter `r` passed through the existing
`ep` common block (as Kartvelishvili's `alpha` is):

```fortran
      elseif(ifrag.eq.4) then
c BCFY vector S-wave heavy-light fragmentation model.
         r=ep
         den=1.d0-(1.d0-r)*z
         poly=2.d0-2.d0*(3.d0-2.d0*r)*z
     #        +3.d0*(3.d0-2.d0*r+4.d0*r**2)*z**2
     #        -2.d0*(1.d0-r)*(4.d0-r+2.d0*r**2)*z**3
     #        +(1.d0-r)**2*(3.d0-2.d0*r+2.d0*r**2)*z**4
         fragfun=xnorm*3.d0*r*z*(1.d0-z)**2*poly/den**6
      elseif(ifrag.eq.5.or.ifrag.eq.8) then
c BCFY pseudoscalar S-wave heavy-light fragmentation model.
         r=ep
         den=1.d0-(1.d0-r)*z
         poly=6.d0-18.d0*(1.d0-2.d0*r)*z
     #        +(21.d0-74.d0*r+68.d0*r**2)*z**2
     #        -2.d0*(1.d0-r)*(6.d0-19.d0*r+18.d0*r**2)*z**3
     #        +3.d0*(1.d0-r)**2*(1.d0-2.d0*r+2.d0*r**2)*z**4
         fragfun=xnorm*r*z*(1.d0-z)**2*poly/den**6
```

Option 8 is an alias of 5 for drivers that label pseudoscalar species
separately. Option 4 is used for the D* feeddown stream of the charm grid,
option 5 for its direct D0 stream; see `../README.md`.

## Terms

FONLL is the work of M. Cacciari, S. Frixione and P. Nason and is
distributed by its authors; only the modifications in these diffs are
part of this repository. Obtain FONLL from the official distribution and
keep its attribution and citation requirements intact.
