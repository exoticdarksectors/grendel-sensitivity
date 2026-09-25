"""Analysis-level constants shared by every benchmark model."""

# HL-LHC integrated luminosity
L_INT_FB = 3000.0           # fb-1
L_INT_PB = L_INT_FB * 1e3   # pb-1

# Exclusion threshold on the expected signal yield. N_signal >= 3 is the
# zero-background 95% CL Poisson upper limit for an observed zero
# (e^-3 ~ 0.05). With no background the median expected observation is
# zero, so this is also the median-expected limit. The selection cuts are not
# here: they live in ``grendel.reco.acceptance``, single-sourced with the
# detector reconstruction.
N_THRESHOLD = 3.0

# Default coupling-scan range (log10 space) for the mixing-angle models.
LOG_U2_MIN = -12.0
LOG_U2_MAX = -1.0
N_U2_POINTS = 200

# Interaction point for the ray-cast: the CMS IP5 origin.
CMS_ORIGIN = (0.0, 0.0, 0.0)

# Flavors known to the HNL analysis and plotting code.
FLAVORS = ["Ue", "Umu", "Utau"]
DEFAULT_FLAVORS = ["Umu"]
