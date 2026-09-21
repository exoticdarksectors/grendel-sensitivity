from grendel.models.scalar import production
from grendel.models.scalar.spec import reconstruction_seed as _reconstruction_seed


def test_reconstruction_seed_is_stable_for_canonical_subsets():
    for index in (0, len(production.MASS_GRID) // 2, len(production.MASS_GRID) - 1):
        mass = production.MASS_GRID[index]
        assert _reconstruction_seed(mass) == 1000 + index
