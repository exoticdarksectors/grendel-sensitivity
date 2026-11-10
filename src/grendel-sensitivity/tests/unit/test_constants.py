"""The run conditions and the exclusion criterion every model shares."""
import numpy as np

from grendel import constants


def test_luminosity_and_threshold():
    assert constants.L_INT_FB == 3000.0
    assert constants.L_INT_PB == 1000.0 * constants.L_INT_FB
    assert constants.N_THRESHOLD == 3.0


def test_coupling_grid_is_log_spaced_over_eleven_decades():
    grid = np.logspace(constants.LOG_U2_MIN, constants.LOG_U2_MAX, constants.N_U2_POINTS)
    assert (constants.LOG_U2_MIN, constants.LOG_U2_MAX, constants.N_U2_POINTS) == (-12.0, -1.0, 200)
    assert np.allclose(np.diff(np.log10(grid)), np.diff(np.log10(grid))[0])


def test_flavours_and_origin():
    assert tuple(constants.FLAVORS) == ("Ue", "Umu", "Utau")
    assert set(constants.DEFAULT_FLAVORS) <= set(constants.FLAVORS)
    assert tuple(constants.CMS_ORIGIN) == (0.0, 0.0, 0.0)
