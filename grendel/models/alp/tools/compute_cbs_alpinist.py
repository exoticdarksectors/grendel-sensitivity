"""Evaluate the universal-fermion ALP FCNC coefficient c_bs (BC10, Lambda=1 TeV)
with the ALPINIST implementation of the GKOZ (arXiv:2310.03524) RG equations.

This evaluates ALPINIST's ``cbs_t``; the convention-mapped value used as
``model.CBS_EFF`` is printed as ``cbs_t/2``. It also prints the low-scale
universal charged-lepton coefficient as the diagnostic ``cel``. The current
decay model does not consume a separate ``CLL_RG`` constant. Run:

    python alp_fermion/tools/compute_cbs_alpinist.py

Everything below is a verbatim port of ALP_rescale/general/{constants,functions}.py
and ALP_rescale/alp/effective_coupling.py (class above_EW) from
github.com/jjerhot/ALPINIST (BSD 3-Clause; Copyright (c) 2022 Jan Jerhot,
Babette Dobrich, Fatih Ertas, Felix Kahlhoefer, Tommaso Spadaro), trimmed to
what the b->s and lepton coefficients need.  Pinned source commit: see
../data/alpinist/COMMIT_SHA.txt.  The universal-fermion scenario (C_f = 0.5,
C_F = -0.5, C_tt = 1) follows ALP_rescale/general/load_data.py.
"""
import numpy as np
from scipy.integrate import quad

# --- constants (ALP_rescale/general/constants.py) ---
m_Z = 91.1876
m_W = 80.379
theta_w = np.arccos(m_W / m_Z)
sin_w2 = 1 - (m_W / m_Z) ** 2
m_u, m_d, m_s, m_c, m_b, m_t = 2.16e-3, 4.7e-3, 93.5e-3, 1.273, 4.183, 172.57
m_q = [m_u, m_d, m_s, m_c, m_b, m_t]
m_qt = [m_u, m_c, m_t]
m_qb = [m_d, m_s, m_b]
m_el, m_mu, m_tau = 5.1099895e-4, 0.1056583755, 1.77693
m_l = [m_el, m_mu, m_tau]
V_tb, V_ts, V_td = 0.99917, 0.03978, 0.0085
V_ub, V_cb = 0.00361, 0.0405

# --- functions (ALP_rescale/general/functions.py) ---
def alpha_s_perturbative(q, order=0, nf=-1, Lambda_QCD=0.34):
    beta_0 = 11. - 2 / 3. * nf
    log_mL = np.log((q / Lambda_QCD) ** 2)
    alpha_LO = 4 * np.pi / (beta_0 * log_mL)
    if order == 0:
        return alpha_LO
    log_log_mL = np.log(log_mL)
    beta_1 = 102. - 38 / 3 * nf
    if order == 1:
        return alpha_LO * (1 - beta_1 / beta_0 ** 2 * log_log_mL / log_mL)
    beta_2 = 2857 / 2 - 5033 / 18 * nf + 325 / 54 * nf ** 2
    if order == 2:
        return alpha_s_perturbative(q, 1, nf=nf, Lambda_QCD=Lambda_QCD) + alpha_LO * (
            beta_1 ** 2 / (beta_0 ** 4 * log_mL ** 2) * (log_log_mL ** 2 - log_log_mL - 1 + beta_2 * beta_0 / beta_1 ** 2))
    xi = 1.2026
    beta_3 = (149753 / 6 + 3564 * xi) - (1078361 / 162 + 6508 / 27 * xi) * nf + (50065 / 162 + 6472 / 81 * xi) * nf ** 2 + 1093 / 729 * nf ** 3
    return alpha_s_perturbative(q, 2, nf=nf, Lambda_QCD=Lambda_QCD) + alpha_LO * (
        beta_1 ** 3 / (beta_0 ** 6 * log_mL ** 3) * (-log_log_mL ** 3 + 5 / 2 * log_log_mL ** 2 + 2 * log_log_mL - 1 / 2
                                                     - 3 * beta_2 * beta_0 / beta_1 ** 2 * log_log_mL + beta_3 * beta_0 ** 2 / (2 * beta_1 ** 3)))

def adaptive_nf_LambdaQCD(q):
    nf, Lambda_QCD = 2, 0.31
    for mq in [m_s, m_c, m_b, m_t]:
        if mq < q:
            nf += 1
            Lambda_QCD = Lambda_QCD * (Lambda_QCD / mq) ** (2 / (33 - 2 * nf))
    return nf, Lambda_QCD

def alpha_s(q, order=3, q_suspension=[0, 1], nf=-1, Lambda_QCD=0.34):
    if q < q_suspension[0]:
        return 1
    if q >= q_suspension[1]:
        if nf == -1:
            nf, Lambda_QCD = adaptive_nf_LambdaQCD(q)
        return alpha_s_perturbative(q, order=order, nf=nf, Lambda_QCD=Lambda_QCD)
    if nf == -1:
        nf, Lambda_QCD = adaptive_nf_LambdaQCD(q_suspension[1])
    x1, y1, m1 = q_suspension[0], 1, 0
    x2 = q_suspension[1]
    y2 = alpha_s_perturbative(x2, order=order, nf=nf, Lambda_QCD=Lambda_QCD)
    m2 = (alpha_s_perturbative(1.000001 * x2, order=order, nf=nf, Lambda_QCD=Lambda_QCD) - y2) * 1000000 / x2
    a_fit = (m1 + m2 - 2 * (y2 - y1) / (x2 - x1)) / (x1 - x2) ** 2
    b_fit = (m2 - m1) / (2 * (x2 - x1)) - (3 / 2) * (x1 + x2) * a_fit
    c_fit = m1 - 3 * (x1 ** 2) * a_fit - 2 * x1 * b_fit
    d_fit = y1 - (x1 ** 3) * a_fit - (x1 ** 2) * b_fit - x1 * c_fit
    return a_fit * q ** 3 + b_fit * q ** 2 + c_fit * q + d_fit

def alpha_t(mu, order=3, q_suspension=[0, 1], nf=-1, Lambda_QCD=0.34):
    alpha_s_ref = alpha_s(m_q[5], order, q_suspension, nf, Lambda_QCD)
    alpha_s_mu = alpha_s(mu, order, q_suspension, nf, Lambda_QCD)
    return 1. / (4 * np.pi) * (alpha_s_mu / alpha_s_ref) ** (8 / 7) / (1 + 9 / 2 * (4 * np.pi * alpha_s_ref) * ((alpha_s_mu / alpha_s_ref) ** (1 / 7) - 1))

def alpha_ww(mu):
    alpha_EM_Z = 1. / 127.952
    alpha_ww_ref = alpha_EM_Z / np.sin(theta_w) ** 2
    beta_ww = 19 / 6
    return alpha_ww_ref / (1 - beta_ww / (2 * np.pi) * alpha_ww_ref * np.log(m_Z / mu))

def alpha_bb(mu):
    alpha_EM_Z = 1. / 127.952
    alpha_bb_ref = alpha_EM_Z / np.cos(theta_w) ** 2
    beta_bb = -41 / 6
    return alpha_bb_ref / (1 - beta_bb / (2 * np.pi) * alpha_bb_ref * np.log(m_Z / mu))

def beta_QED(mu):
    beta = 0
    for m in m_qt:
        if mu > m:
            beta += 3 * (2. / 3) ** 2
    for m in m_qb:
        if mu > m:
            beta += 3 * (-1. / 3) ** 2
    for m in m_l:
        if mu > m:
            beta += 1
    return -4. / 3 * beta

def alpha_EM(mu):
    alpha_EM_Z = 1 / 127.952
    return alpha_EM_Z / (1 - beta_QED(mu) / (2 * np.pi) * alpha_EM_Z * np.log(m_Z / mu))

# --- above_EW (ALP_rescale/alp/effective_coupling.py) ---
class above_EW:
    def __init__(self, matching_scale=1000.):
        self._lambda = matching_scale

    def C_GG_tilde(self, C_GG, C_u=[0, 0, 0], C_d=[0, 0, 0], C_Q=[0, 0, 0]):
        return C_GG + 1. / 2 * (np.array(C_u) + np.array(C_d) - 2 * np.array(C_Q)).sum()

    def C_WW_tilde(self, C_WW, C_Q=[0, 0, 0], C_L=[0, 0, 0]):
        return C_WW - 1. / 2 * (3 * np.array(C_Q) + np.array(C_L)).sum()

    def C_BB_tilde(self, C_BB, C_u=[0, 0, 0], C_d=[0, 0, 0], C_Q=[0, 0, 0], C_e=[0, 0, 0], C_L=[0, 0, 0]):
        return C_BB + (4. / 3 * np.array(C_u) + 1. / 3 * np.array(C_d) - 1. / 6 * np.array(C_Q) + np.array(C_e) - 1. / 2 * np.array(C_L)).sum()

    def __one_minus_e18U(self, mu, lam):
        return 9. / 2 * alpha_t(mu) / alpha_s(mu) * (1 - (alpha_s(lam) / alpha_s(mu)) ** (1 / 7))

    def C_GG_tilde_mu(self, mu, C_GG_tilde, C_tt):
        return C_GG_tilde - 2. / 9 * self.__one_minus_e18U(mu, self._lambda) * C_tt

    def C_WW_tilde_mu(self, mu, C_WW_tilde, C_tt):
        return C_WW_tilde - 1. / 6 * self.__one_minus_e18U(mu, self._lambda) * C_tt

    def C_BB_tilde_mu(self, mu, C_BB_tilde, C_tt):
        return C_BB_tilde - 17. / 54 * self.__one_minus_e18U(mu, self._lambda) * C_tt

    def __I_t(self, mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt):
        def integrand_Lambda_mu(mu):
            return 1. / mu * self.__one_minus_e18U(mu_w, mu) \
                * (2. * alpha_s(mu) ** 2 / (np.pi ** 2) * self.C_GG_tilde_mu(mu, C_GG_tilde, C_tt)
                   + 9. * alpha_ww(mu) ** 2 / (16. * np.pi ** 2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt)
                   + 17. * alpha_bb(mu) ** 2 / (48. * np.pi ** 2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return -2. / 3 * self.__one_minus_e18U(mu_w, self._lambda) * C_tt - 2. / 3 * integral_Lambda_mu[0]

    def C_uu_mu(self, mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt, C_uu=None, type=None):
        delta = 1 if type == "t" else 0
        if C_uu is None:
            C_uu = C_tt
        def integrand_Lambda_mu(mu):
            return 1. / mu * (1 - 2. / 3 * (1 + delta / 2) * self.__one_minus_e18U(mu_w, mu)) \
                * (2. * alpha_s(mu) ** 2 / np.pi ** 2 * self.C_GG_tilde_mu(mu, C_GG_tilde, C_tt)
                   + 9. * alpha_ww(mu) ** 2 / (16. * np.pi ** 2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt)
                   + 17. * alpha_bb(mu) ** 2 / (48. * np.pi ** 2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return C_uu - 2. / 3 * (1 + delta / 2) * self.__one_minus_e18U(mu_w, self._lambda) * C_tt + integral_Lambda_mu[0]

    def C_ee_mu(self, mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt, C_ee=None):
        if C_ee is None:
            C_ee = C_tt
        def integrand_Lambda_mu(mu):
            return 1. / mu * (1 - self.__one_minus_e18U(mu_w, mu)) \
                * (9. * alpha_ww(mu) ** 2 / (16. * np.pi ** 2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt)
                   + 15. * alpha_bb(mu) ** 2 / (16. * np.pi ** 2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return C_ee - self.__I_t(mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt) + integral_Lambda_mu[0]

    def C_qq_ij_mu(self, mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt, C_qq_ij=0, i=None, j=None):
        types = {"t": V_tb, "c": V_cb, "u": V_ub, "b": V_tb, "s": V_ts, "d": V_td}
        V_i, V_j = types[i], types[j]
        x_t = m_t ** 2 / m_W ** 2
        return C_qq_ij + V_i * V_j * (
            -1. / 6 * self.__I_t(mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt)
            + alpha_t(mu_w) / (4 * np.pi) * (
                self.C_uu_mu(mu_w, C_GG_tilde, C_WW_tilde, C_BB_tilde, C_tt, type="t")
                * (np.log(mu_w ** 2 / m_t ** 2) / 2 - 1. / 4 - 3. / 2 * (1 - x_t + np.log(x_t)) / (1 - x_t) ** 2)
                - 3 * alpha_EM(mu_w) / (2 * np.pi * sin_w2) * self.C_WW_tilde_mu(mu_w, C_WW_tilde, C_tt)
                * (1 - x_t + x_t * np.log(x_t)) / (1 - x_t) ** 2))


if __name__ == "__main__":
    LAM = 1000.0
    C_f = [0.5, 0.5, 0.5]
    C_F = [-0.5, -0.5, -0.5]
    above = above_EW(LAM)
    cgg_t = above.C_GG_tilde(0, C_u=C_f, C_d=C_f, C_Q=C_F)
    cww_t = above.C_WW_tilde(0, C_Q=C_F, C_L=C_F)
    cbb_t = above.C_BB_tilde(0, C_u=C_f, C_d=C_f, C_Q=C_F, C_e=C_f, C_L=C_F)
    print(f"C_GG_tilde = {cgg_t:.6f}, C_WW_tilde = {cww_t:.6f}, C_BB_tilde = {cbb_t:.6f}")

    cbs_t = above.C_qq_ij_mu(m_t, cgg_t, cww_t, cbb_t, 1, C_qq_ij=0, i="b", j="s")
    print(f"cbs_t (universal fermion, Lambda=1 TeV)         = {cbs_t:.6e}")
    print(f"  -> g_bs = cbs_t/Lambda                        = {cbs_t / LAM:.6e} GeV^-1")

    cel = above.C_ee_mu(m_t, cgg_t, cww_t, cbb_t, 1)
    print(f"cel = cmu = ctau (C_ee_mu at mu_w=m_t)          = {cel:.6f}")

    # leading-log comparison (our current model.h_sb, conventions: amplitude uses h_sb*inv_f)
    V_HIGGS = 246.21965
    M_TOP_OLD = 172.5
    V_TB_VTS = 0.0405
    h_sb_old = 1 / (16 * np.pi ** 2) * (M_TOP_OLD ** 2 / V_HIGGS ** 2) * V_TB_VTS * np.log(LAM ** 2 / M_TOP_OLD ** 2)
    print(f"old leading-log h_sb                            = {h_sb_old:.6e}")
    print(f"ALPINIST equivalent in our convention cbs_t/2   = {cbs_t / 2:.6e}")
    print(f"ratio (new/old) amplitude                        = {cbs_t / 2 / h_sb_old:.4f}")
    print(f"ratio (new/old) rate                             = {(cbs_t / 2 / h_sb_old) ** 2:.4f}")
