"""
hnl/analysis/constants.py

Analysis-level constants for the GRENDEL HNL sensitivity scan.

DEFAULT_FLAVORS is the muon-mixing scenario. ANALYSIS_MASS_MAX is 10 GeV so
the committed electroweak production tail is included.
"""

# HL-LHC integrated luminosity
L_INT_FB = 3000.0           # fb-1
L_INT_PB = L_INT_FB * 1e3   # pb-1 = 3e6 pb-1

# Exclusion threshold on the expected signal yield. N_signal >= 3 is the
# zero-background 95% CL Poisson upper limit for an observed zero
# (e^-3 ~ 0.05), the community convention for HNL money plots. (The
# median-expected zero-background limit is ~2.44; 3.0 is the slightly
# conservative observed-zero number.) The live selection cuts are NOT here:
# they live in analysis/decay_reco_acceptance.py, single-sourced with the
# higgs/ GRENDEL reconstruction.
N_THRESHOLD = 3.0

# U^2 scan range (log10 space)
LOG_U2_MIN = -12.0
LOG_U2_MAX = -1.0
N_U2_POINTS = 200

# CMS IP5 origin for ray-casting
CMS_ORIGIN = (0.0, 0.0, 0.0)

# Full analysis range represented by the committed mass grid.
ANALYSIS_MASS_MAX = 10.0  # GeV

# Flavors known to the analysis and plotting code.
FLAVORS = ["Ue", "Umu", "Utau"]

# Production focus: muon-only scenario (HNL pattern 010).
DEFAULT_FLAVORS = ["Umu"]
