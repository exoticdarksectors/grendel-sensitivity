#!/usr/bin/env bash
# 01_fetch_sources.sh -- fetch the Mathematica-side sources at their pins (no Wolfram kernel needed).
#
#   source env/activate.sh && bash pipelines/senscalc/01_fetch_sources.sh
#
#  * SensCalc (github.com/maksymovchynnikov/SensCalc), shallow fetches of two commits into one git repo:
#      third_party/SensCalc         eb2c456e8c9b35e77d212d18813306da195a13cb  (2026-10-17; FASER2-FPF / SHiP BC10 curves)
#      third_party/SensCalc-v1.4    4f7210b6e566f5510e2c7269adec3744a6369880  (tag v1.4, 2026-10-11; senscalc_2501 export,
#                                   a git worktree of the same repo, because export_senscalc_2501.py checks `git rev-parse HEAD`)
#  * FeynCalc 10.2.1 (tag Release-10_2_1, the version of the reference runs), unpacked to env/wolfram/Applications/FeynCalc.
#    SensCalc's notebooks 1 and 2 do Get["FeynCalc`"]; our drivers put env/wolfram/Applications on $Path.
#
# Both SensCalc checkouts stay pristine: every run works in a `git archive` copy under work/senscalc/.
# Validation (logs/senscalc/01_fetch_sources.log):
#  * every file of the eb2c456 tree has the git blob hash of the member of the GitHub tarball the campaign used
#    (campaigns/bc10/faser2_senscalc_20261025/src/SensCalc-eb2c456….tar.gz, sha256 de29db62…d705, in the optional
#    reference tree $GRENDEL_REF; skipped without it);
#  * the v1.4 MX inputs match the hashes pinned in export_senscalc_2501*.py (--check-only, no kernel).
# Resume-safe: existing checkouts at the right commit are kept.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"
source "$REPRO/pipelines/common/ref.sh"

URL=https://github.com/maksymovchynnikov/SensCalc
EB=eb2c456e8c9b35e77d212d18813306da195a13cb
V14=4f7210b6e566f5510e2c7269adec3744a6369880
SC="$REPRO/third_party/SensCalc"
SC14="$REPRO/third_party/SensCalc-v1.4"
CAMPAIGN_TGZ=$REF/grendel/campaigns/bc10/faser2_senscalc_20261025/src/SensCalc-eb2c456e8c9b35e77d212d18813306da195a13cb.tar.gz
FC_TAG=Release-10_2_1
FC_COMMIT=782aa16e71          # prefix of the tag commit (GitHub API, 2026-11-06)
FC_DIR="$REPRO/env/wolfram/Applications/FeynCalc"
DL="$REPRO/work/senscalc/downloads"
mkdir -p "$DL" "$REPRO/env/wolfram/Applications"

# ---------------- SensCalc ----------------
if [ ! -d "$SC/.git" ]; then
  git init -q "$SC"
  git -C "$SC" remote add origin "$URL"
fi
for c in $EB $V14; do
  git -C "$SC" cat-file -e "$c^{commit}" 2>/dev/null || git -C "$SC" fetch --depth 1 origin "$c"
done
[ "$(git -C "$SC" rev-parse HEAD 2>/dev/null || true)" = "$EB" ] || git -C "$SC" checkout -q --detach "$EB"
if [ ! -e "$SC14/.git" ]; then
  git -C "$SC" worktree add -q --detach "$SC14" "$V14"
fi
[ "$(git -C "$SC14" rev-parse HEAD)" = "$V14" ] || git -C "$SC14" checkout -q --detach "$V14"
echo "SensCalc eb2c456: $(git -C "$SC" rev-parse HEAD)  $(du -sh --exclude=.git "$SC" | cut -f1)"
echo "SensCalc v1.4   : $(git -C "$SC14" rev-parse HEAD)"
git -C "$SC" status --porcelain | head -5

# same tree as the GitHub tarball the campaign extracted?  (git blob hashes of every member; read-only on the reference)
if have_ref "SensCalc eb2c456 checkout vs the campaign's GitHub tarball" "$CAMPAIGN_TGZ"; then
  git -C "$SC" ls-tree -r "$EB" > "$REPRO/work/senscalc/senscalc_eb2c456_lstree.txt"
  python3 -I "$REPRO/pipelines/senscalc/py/compare_tarball_tree.py" "$CAMPAIGN_TGZ" "$REPRO/work/senscalc/senscalc_eb2c456_lstree.txt"
fi

# pinned v1.4 inputs of the senscalc_2501 exporters (hash check only, no kernel)
( cd "$REPRO/src/grendel-sensitivity" && python -m grendel.models.alp.tools.export_senscalc_2501 "$SC14" --check-only \
    && python -m grendel.models.alp.tools.export_senscalc_2501_production "$SC14" --check-only )

# ---------------- FeynCalc ----------------
if [ ! -f "$FC_DIR/FeynCalc.m" ]; then
  tgz="$DL/feyncalc-$FC_TAG.tar.gz"
  [ -s "$tgz" ] || curl -fsSL --retry 3 -o "$tgz" "https://github.com/FeynCalc/feyncalc/archive/refs/tags/$FC_TAG.tar.gz"
  sha256sum "$tgz"
  tmp=$(mktemp -d "$DL/fc.XXXX")
  tar -xzf "$tgz" -C "$tmp"
  rm -rf "$FC_DIR"; mv "$tmp/feyncalc-$FC_TAG/FeynCalc" "$FC_DIR"; rm -rf "$tmp"
fi
grep -m1 -E 'FeynCalcVersion|\$FeynCalcVersion *=' "$FC_DIR"/FeynCalc.m || true
echo "FeynCalc in $FC_DIR"
