"""The HNL mass grids of the published v2.3 table (curves/data/grendel_v2_3/hnl, 241 masses per flavour).

COARSE   the 104 masses per flavour of the as-built 1x production scanned on lxplus (v2_3_20261016/lxplus/
         points_hnl.txt): the repo MASS_GRID at 0.2-1.975 and 4.4-10 GeV.
FINE     the 137 masses per flavour of hnl_fine_condor_20261014/masses.txt: 2.0-3.575 GeV in 25 MeV,
         3.6-3.9 in 5 MeV, 3.925-4.2 in 25 MeV.
JUNE116  the 116-mass grid of the June 2026 v1 production (llpatcolliders 06e5a18 hnl/config_mass_grid.py): the
         repo MASS_GRID without the 3.62-3.70 closure refinement. Channels generated over a mass list share one
         random stream, so the v1 channels must be generated over exactly this list.
"""
from grendel.models.hnl.mass_grid import MASS_GRID

FLAVORS = ["Ue", "Umu", "Utau"]
COARSE = [m for m in MASS_GRID if m < 2.0 - 1e-9 or m > 4.2 + 1e-9]
FINE = sorted({round(2.0 + 0.025 * i, 3) for i in range(64)} | {round(3.6 + 0.005 * i, 3) for i in range(61)}
              | {round(3.925 + 0.025 * i, 3) for i in range(12)})
JUNE116 = [m for m in MASS_GRID if m not in (3.62, 3.63, 3.64, 3.65, 3.675, 3.685, 3.70)]
SEEDS_FINE = {"Ue": 7000, "Umu": 7100, "Utau": 7200}
# MadGraph seed offsets (SubProcesses/randinit) of the published W/Z samples
RANDINIT_COARSE = 21                       # June v1: fresh process directory (MG5 template r=21) at every mass
RANDINIT_FINE = {"Ue": 27, "Umu": 33, "Utau": 27}
RANDINIT_FINE_EXCEPTIONS = {("Umu", 3.845): 34}   # lxplus cluster 16488236 (three failures with r=33)


def fine_randinit(flavor: str, mass: float) -> int:
    return RANDINIT_FINE_EXCEPTIONS.get((flavor, round(mass, 3)), RANDINIT_FINE[flavor])


def label(m: float) -> str:
    return f"{m:.3f}".replace(".", "p")


assert len(COARSE) == 104 and len(FINE) == 137 and len(JUNE116) == 116
assert not set(COARSE) & set(FINE)
