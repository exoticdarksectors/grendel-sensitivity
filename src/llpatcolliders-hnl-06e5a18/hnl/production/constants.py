"""Physical constants for HNL production at LHC 14 TeV."""

from particle import Particle

M_B0 = Particle.from_pdgid(511).mass * 1e-3      # B0
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3    # B+
M_BS = Particle.from_pdgid(531).mass * 1e-3       # Bs
M_BC = Particle.from_pdgid(541).mass * 1e-3       # Bc+
M_D0 = Particle.from_pdgid(421).mass * 1e-3       # D0
M_DPLUS = Particle.from_pdgid(411).mass * 1e-3    # D+
M_DS = Particle.from_pdgid(431).mass * 1e-3       # Ds+
M_TAU = Particle.from_pdgid(15).mass * 1e-3       # tau
M_ELECTRON = Particle.from_pdgid(11).mass * 1e-3  # electron
M_MUON = Particle.from_pdgid(13).mass * 1e-3      # muon
M_PION = Particle.from_pdgid(211).mass * 1e-3     # pi+
M_KAON = Particle.from_pdgid(321).mass * 1e-3     # K+
M_DSTAR0 = Particle.from_pdgid(423).mass * 1e-3   # D*0   (B -> D* tau nu daughter)
M_DSTARP = Particle.from_pdgid(413).mass * 1e-3   # D*+   (B -> D* tau nu daughter)
M_DSSTAR = Particle.from_pdgid(433).mass * 1e-3   # Ds*+  (Bs -> Ds* tau nu daughter)
M_LAMBDA_B = Particle.from_pdgid(5122).mass * 1e-3  # Lambda_b0
M_LAMBDA_C = Particle.from_pdgid(4122).mass * 1e-3  # Lambda_c+

# Inclusive pp -> Bc+ + Bc-; already both charges, so no downstream factor 2.
# This is a model input, not a measured inclusive 14 TeV cross section. Replace
# it together with the borrowed bottom shape when a dedicated Bc grid is
# supplied (`REMAINING_WORK.md`).
SIGMA_BC_PB = 0.9e6

# QCD K-factor applied to the LO MadGraph electroweak rates. The single
# inclusive value below is a documented approximation; a process- and
# mass-differential NLO/LO determination requires dedicated NLO MadGraph runs
# (see `REMAINING_WORK.md`). The per-process table is keyed so those
# inputs can be dropped in without touching the drivers; all entries currently
# hold the inclusive default so behaviour is unchanged until real values land.
#
# Only keys that the drivers can actually apply are listed. The W/Z -> l N
# sample mixes W -> l N (dominant) and Z -> nu N in one LHE/CSV with no per-row
# process tag, so it is scaled by the single "W" key; a separate "Z" value
# would have no effect without carrying the LHE mother PDG into the HNL rows (a
# future refinement, see `REMAINING_WORK.md`). The prompt-tau pool is
# separable by parent origin, so "W" and "DY" are applied per row there.
K_FACTOR_EW = 1.3
K_FACTOR_EW_BY_PROCESS = {
    "W": K_FACTOR_EW,   # p p -> W -> l N  and  W -> tau nu (incl. the unsplit Z -> nu N)
    "DY": K_FACTOR_EW,  # p p -> gamma*/Z -> tau tau (Drell-Yan tau pool)
}

# Approximate K+ + K- flux: no current data/generator validation is attached.
# Replace with a measured or Pythia-derived (pT, y) spectrum when available.
SIGMA_KAON_PB = 3.0e11  # ~80 mb × ~4 K± per inelastic event (approximate)

KAON_TSALLIS_T = 0.17   # GeV
KAON_TSALLIS_N = 7.0
KAON_PT_MAX = 5.0       # GeV, sampling ceiling (kaons of interest are soft)
KAON_RAPIDITY_SIGMA = 2.5
KAON_E_MAX = 7000.0  # GeV, half of sqrt(s) at LHC 14 TeV

# Species fractions are acceptance-specific inputs, applied over the full FONLL
# grid as constants. See FRAGMENTATION_POLICY and `REMAINING_WORK.md`.
FRAG_B = {
    521: 0.36205648081100655,  # B+/B-
    511: 0.36205648081100655,  # B0/B0bar
    531: 0.08834178131788560,  # Bs
}

FRAG_C = {
    421: 0.382,   # D0/D0bar, includes strong D* feeddown convention
    411: 0.191,   # D+/D-
    431: 0.061,   # Ds+/Ds-
}

OMITTED_FRAG_B = {
    "b_baryons": 0.18754525706010140,
}

# The closure remainder of the LHCb fu, fd, fs, and Lambda_b ratios is
# represented as Lambda_b-like production. It is not a direct measurement of a
# full b-baryon fraction: Xi_b/Omega_b are neither separated nor validated as
# Lambda_b equivalents. The pT-y shape reuses the bottom FONLL grid.
FRAG_LAMBDA_B = OMITTED_FRAG_B["b_baryons"]

OMITTED_FRAG_C = {
    "Lambda_c+": 0.168,
    "Xi_c0": 0.099,
    "Xi_c+": 0.096,
    "J/psi": 0.0037,
}

FRAGMENTATION_POLICY = {
    "bottom": {
        "source": "LHCb Phys. Rev. D100 (2019) 031102, arXiv:1902.06794",
        "scope": "pp 13 TeV, 2<eta<5, 4<pT<25 GeV; extrapolated over the FONLL grid",
        "derivation": "fu=fd; closure using fs/(fu+fd)=0.122 and fLambda_b/(fu+fd)=0.259",
        "meson_fraction_sum": sum(FRAG_B.values()),
        "omitted_fraction_sum": sum(OMITTED_FRAG_B.values()),
    },
    "charm": {
        "source": "ALICE charm fragmentation fractions in pp 13 TeV, arXiv:2308.04877",
        "scope": "pp 13 TeV, |y|<0.5; extrapolated over the FONLL grid",
        "derivation": "D* feeddown counted in D0/D+ rather than as an independent weak parent",
        "meson_fraction_sum": sum(FRAG_C.values()),
        "omitted_fraction_sum": sum(OMITTED_FRAG_C.values()),
    },
}

MESON_MASSES = {
    521: M_BPLUS, 511: M_B0, 531: M_BS, 541: M_BC,
    421: M_D0, 411: M_DPLUS, 431: M_DS,
}

LEPTON_MASSES = {
    'Ue': M_ELECTRON,
    'Umu': M_MUON,
    'Utau': M_TAU,
}

FLAVOR_TO_LEPTON_PDG = {
    'Ue': 11,
    'Umu': 13,
    'Utau': 15,
}

FLAVOR_TO_MG5 = {
    'Ue': 'electron',
    'Umu': 'muon',
    'Utau': 'tau',
}

QUARK_MESON_MAP = {
    'bottom': [(521, FRAG_B[521]), (511, FRAG_B[511]), (531, FRAG_B[531])],
    'charm':  [(421, FRAG_C[421]), (411, FRAG_C[411]), (431, FRAG_C[431])],
}

# Parent -> tau X normalization inputs for the induced-tau pool. These are
# central values only; uncertainties are not yet propagated.
#
# PDG 2025 update:
BR_DS_TAUNU        = 5.39e-2   # Ds+ -> tau+ nu
BR_DPLUS_TAUNU     = 1.20e-3   # D+  -> tau+ nu
BR_BPLUS_TAUNU     = 1.09e-4   # B+  -> tau+ nu
BR_BP_D0_TAUNU     = 0.77e-2   # B+  -> D0bar  tau+ nu
BR_BP_DSTAR0_TAUNU = 1.88e-2   # B+  -> D*0bar tau+ nu
BR_B0_DP_TAUNU     = 0.98e-2   # B0  -> D-     tau+ nu
BR_B0_DSTARP_TAUNU = 1.48e-2   # B0  -> D*-    tau+ nu

# SM lattice-derived central values. The Bs values combine HPQCD
# R(Ds), R(Ds*) and light-lepton normalization results (arXiv:1906.00701,
# arXiv:2105.11433); neither mode has a direct branching-fraction measurement.
BR_BS_DS_TAUNU     = 0.71e-2   # Bs -> Ds-  tau+ nu
BR_BS_DSSTAR_TAUNU = 1.33e-2   # Bs -> Ds*- tau+ nu

# Theory estimate from the leptonic-decay formula (f_Bc and Vcb dependent);
# there is no direct Bc -> tau nu measurement.
BR_BC_TAUNU        = 2.3e-2

# Direct LHCb measurement, arXiv:2201.03497; its third uncertainty comes from
# the external normalization branching fraction.
BR_LB_LC_TAUNU     = 1.50e-2
