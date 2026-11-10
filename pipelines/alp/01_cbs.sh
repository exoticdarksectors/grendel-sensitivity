#!/usr/bin/env bash
# Step 01 -- the b -> s a coupling C_bs behind the BC10 production weights.
#   * grendel/models/alp/tools/compute_cbs_alpinist.py: the ALPINIST port of the GKOZ RG equations (cross-check;
#     the v1/v2 value 3.518383e-4);
#   * model.CBS_EFF = C_bs f_GKOZ / (m_b - m_s) from GKOZ (arXiv:2310.03524) Table 1, |C_bs| = 1.8e-9 at
#     f_GKOZ = 1 PeV: 4.405e-4, used since v2.2 (CHANGELOG 2.2, curves/data/grendel_v2_2/README.md);
#   * the v2.2 weight factor (4.405/3.518)^2 = 1.567236 and BR(B+ -> K+ a) at f_GKOZ = 1 PeV (GKOZ Fig. 4: 1.5e-7).
# Output: work/alp/01_cbs/cbs.txt. Seconds; no slot needed.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
OUT="$W/01_cbs"; mkdir -p "$OUT"
{
  echo "# compute_cbs_alpinist.py (ALPINIST port)"
  python "$SRC/grendel/models/alp/tools/compute_cbs_alpinist.py"
  echo
  echo "# grendel.models.alp.model"
  python - <<'EOF'
from grendel.models.alp import model as m
alpinist = 3.518383e-4
print(f"CBS_EFF = C_bs f_GKOZ/(m_b-m_s) = {m.C_BS_GKOZ_1PEV:g}*1e6/({m.QUARK_MASSES['b']}-{m.QUARK_MASSES['s']}) = {m.CBS_EFF:.6e}")
print(f"(CBS_EFF/ALPINIST)^2 = {(m.CBS_EFF / alpinist) ** 2:.6f}   (v2.2 README: x1.567236)")
inv_f = 2.0 / 1.0e6          # 1/f = 2/f_GKOZ, f_GKOZ = 1 PeV
for ma in (0.3, 1.0, 2.0):
    print(f"BR(B+ -> K+ a) at f_GKOZ = 1 PeV, m_a = {ma} GeV: {m.br_B_to_Ka(ma, inv_f, parent='B+', kaon='K'):.3e}")
assert abs(m.CBS_EFF - 4.405e-4) / 4.405e-4 < 1e-3
assert abs((m.CBS_EFF / alpinist) ** 2 - 1.567236) < 1e-6
EOF
} | tee "$OUT/cbs.txt"
grep -q "cbs_t (universal fermion, Lambda=1 TeV)         = 3.518383e-04" "$OUT/cbs.txt"
echo "01_cbs OK"
