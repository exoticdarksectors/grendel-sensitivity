"""Physical constants for HNL production at LHC 14 TeV."""

from particle import Particle

M_B0 = Particle.from_pdgid(511).mass * 1e-3
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3
M_BS = Particle.from_pdgid(531).mass * 1e-3
M_BC = Particle.from_pdgid(541).mass * 1e-3
M_D0 = Particle.from_pdgid(421).mass * 1e-3
M_DPLUS = Particle.from_pdgid(411).mass * 1e-3
M_DS = Particle.from_pdgid(431).mass * 1e-3
M_TAU = Particle.from_pdgid(15).mass * 1e-3
M_ELECTRON = Particle.from_pdgid(11).mass * 1e-3
M_MUON = Particle.from_pdgid(13).mass * 1e-3
M_PION = Particle.from_pdgid(211).mass * 1e-3
M_KAON = Particle.from_pdgid(321).mass * 1e-3
M_DSTAR0 = Particle.from_pdgid(423).mass * 1e-3
M_DSTARP = Particle.from_pdgid(413).mass * 1e-3
M_DSSTAR = Particle.from_pdgid(433).mass * 1e-3
M_LAMBDA_B = Particle.from_pdgid(5122).mass * 1e-3
M_LAMBDA_C = Particle.from_pdgid(4122).mass * 1e-3

SIGMA_BC_PB = 0.9e6
SIGMA_BC_REL_UNCERT = 0.40

K_FACTOR_EW = 1.3
K_FACTOR_EW_BY_PROCESS = {
    "W": K_FACTOR_EW,
    "DY": K_FACTOR_EW,
}

SIGMA_KAON_PB = 6.535e11
SIGMA_KAON_PB_TSALLIS = 3.0e11

KAON_D_ESC = 1.5
KAON_D_ESC_RANGE = (1.0, 3.0)

KAON_TSALLIS_T = 0.17
KAON_TSALLIS_N = 7.0
KAON_PT_MAX = 5.0
KAON_RAPIDITY_SIGMA = 2.5
KAON_E_MAX = 7000.0

FRAG_B = {
    521: 0.36205648081100655,
    511: 0.36205648081100655,
    531: 0.08834178131788560,
}

FRAG_C = {
    421: 0.382,
    411: 0.191,
    431: 0.061,
}

OMITTED_FRAG_B = {
    "b_baryons": 0.18754525706010140,
}

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

BR_DS_TAUNU        = 5.39e-2
BR_DPLUS_TAUNU     = 1.20e-3
BR_BPLUS_TAUNU     = 1.09e-4
BR_BP_D0_TAUNU     = 0.77e-2
BR_BP_DSTAR0_TAUNU = 1.88e-2
BR_B0_DP_TAUNU     = 0.98e-2
BR_B0_DSTARP_TAUNU = 1.48e-2

BR_BS_DS_TAUNU     = 0.71e-2
BR_BS_DSSTAR_TAUNU = 1.33e-2

BR_BC_TAUNU        = 2.3e-2

BR_LB_LC_TAUNU     = 1.50e-2
