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
# Anchored to LHCb f(Bc)/(f(B-)+f(B0bar)) ~= 3.7e-3 at 13 TeV (arXiv:1910.13404,
# PRD 100 112006): with sigma(B+) ~ 2e8 pb this gives sigma(Bc) ~ 0.8e6 pb, so
# the value below is LHCb-consistent (~16% high). The pT-y shape is still
# borrowed from bottom (a dedicated Bc grid is item 9 in `REMAINING_WORK.md`).
SIGMA_BC_PB = 0.9e6
# Bc production-normalization uncertainty (the B4 nuisance). Conservative for a
# projection: ~25% from BR(Bc->J/psi mu nu) theory folded into the LHCb
# fraction + its pT dependence, in quadrature with ~20% from the borrowed
# bottom shape (Bc is heavier -> harder spectrum -> different acceptance).
SIGMA_BC_REL_UNCERT = 0.40

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

# Charged-kaon flux normalization sigma_inel * <n_K+->. The DEFAULT is the
# Pythia 8.315 SoftQCD:inelastic value at 14 TeV (sigma_inel = 78.93 mb,
# <n_K+-> = 8.28 per inelastic event at the pinned seed 42), measured with
# production/decay_engine/kaon_softqcd.cc and stored in
# production/data/kaon_softqcd_spectrum.npz (regenerate/verify with
# grendel.models.hnl.production.make_kaon_spectrum; requires Pythia
# 8.315 (third_party/fetch.py pythia8315) -- 8.317 shifts <n_K+-> by ~0.5%). The old 3.0e11 stub assumed ~4 K+-
# per inelastic event and is retained for the legacy `--spectrum tsallis` path.
SIGMA_KAON_PB = 6.535e11        # Pythia SoftQCD, seed 42 (supersedes the 3.0e11 stub)
SIGMA_KAON_PB_TSALLIS = 3.0e11  # legacy stub normalization (~80 mb x ~4 K+-)

# Charged-kaon transport: a charged kaon (ctau = 3.712 m) must decay before it is
# absorbed in dense material to produce an escaping HNL. The survival probability
# is P(decay within d_esc) = 1 - exp(-d_esc / (beta*gamma * ctau)), applied as a
# per-kaon weight (the HNL is cast from the IP, not the displaced kaon-decay point;
# decay_engine/transport_control.py measures a displaced/IP acceptance ratio of ~1.0
# within a few % on the sensitivity-relevant long-lifetime plateau, so no per-origin
# change is applied). KAON_D_ESC is the escape path length
# before dense material -- a proxy for the CMS material budget (calorimeter front
# ~1.3 m; the tracker is largely transparent to a decaying kaon). The default and
# the declared [1, 3] m range (KAON_D_ESC_RANGE) are the dominant kaon-sector
# uncertainty; the range is published as the BC6/BC7 transport band
# data/published/bundle/kaon_desc_band.csv (2026-07-18), and d_esc itself remains a
# proxy pending a real material map. Set KAON_D_ESC = None (or --no-transport)
# for the legacy prompt-at-IP behaviour.
KAON_D_ESC = 1.5        # m, central escape distance before dense material
KAON_D_ESC_RANGE = (1.0, 3.0)   # m, published transport-band endpoints

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
