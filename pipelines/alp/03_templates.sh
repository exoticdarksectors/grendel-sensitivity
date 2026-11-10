#!/usr/bin/env bash
# Step 03 -- exHad ALP-fermion rest-frame decay templates (bc10_v24_20261019/make_templates.sh):
#   python -m grendel.models.alp.templates_exhad --n-templates 20000 --seed 1234 --workers 6 --sub-request 250 \
#          --max-retries 10 --exhad-root <exHad e41b07b> --out ... --mass <the 153 masses of points_all.txt>
# run under exHad's venv python (third_party/exHad/.venv, Pythia 8.317) with the grendel code on PYTHONPATH,
# as the original. 151 bundles (0.54 and 0.96 GeV are inside the eta / eta' windows and skipped).
#
#   03_templates.sh [K] [mass ...]     K = parallel single-slot jobs (default 4); masses default to points_all.txt
#
# Each mass has its own seed (sha256 of the base seed and the mass label) and the events of a sub-request depend
# only on its seed and size (sub-request 250 < chunk 512), so a bundle does not depend on the number of exHad
# workers, on the other masses or on the run order (checked: 3 workers here reproduce the 6-worker macOS decays).
# Instead of one K-worker process (a K-slot request starves in the shared slot pool) the masses are dealt
# round-robin to K one-worker jobs of one slot each, each with its own folder work/alp/templates_jobs/job<j>; the
# bundles are then moved into templates_exhad, the jobs' manifest records (with their retry logs) are merged into
# its MANIFEST.json (merge_manifests.py), and a last driver pass over all masses (--skip-existing: nothing left to
# generate) writes the complete manifest.
# Resume: --skip-existing and the per-job folders; bundles that do not open (a job killed while writing) are
# removed first. Output: work/alp/templates_exhad/templates_*.npz + MANIFEST.json. ~0.3 GB RAM per job.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
K="${1:-4}"; shift || true
MASSES="${*:-$(tr '\n' ' ' < "$POINTS")}"
OUT="$W/templates_exhad"; JOBS="$W/templates_jobs"; mkdir -p "$OUT" "$JOBS"

# drop unreadable bundles so that --skip-existing never keeps a truncated file
"$GRENDEL_EXHAD_PYTHON" "$ALP/merge_manifests.py" check "$OUT" "$JOBS"

gen() {   # gen <tag> <out dir> <masses...>
  local tag=$1 out=$2; shift 2
  PYTHONPATH="$SRC" $SLOT -n 1 --label "alp:03_templates_$tag" -- \
    "$GRENDEL_EXHAD_PYTHON" -m grendel.models.alp.templates_exhad \
      --n-templates "$N_TEMPLATES" --seed "$TPL_SEED" --workers 1 --sub-request "$SUB_REQUEST" \
      --max-retries "$MAX_RETRIES" --exhad-root "$GRENDEL_THIRD_PARTY_DIR/exHad" --out "$out" --skip-existing \
      --mass "$@" >> "$L/03_templates_$tag.log" 2>&1
}

echo "$(date '+%F %T') 03_templates start K=$K masses: $MASSES" >> "$L/03_templates.log"
# masses still to do (neither in templates_exhad nor in a job folder), dealt round-robin to K jobs
todo=()
for m in $MASSES; do
  b="templates_$(printf '%.3f' "$m" | tr . p).npz"
  [ -f "$OUT/$b" ] || ls "$JOBS"/job*/"$b" > /dev/null 2>&1 || todo+=("$m")
done
pids=()
for ((j = 0; j < K; j++)); do
  sub=(); for ((i = j; i < ${#todo[@]}; i += K)); do sub+=("${todo[i]}"); done
  [ ${#sub[@]} -gt 0 ] || continue
  gen "job$j" "$JOBS/job$j" "${sub[@]}" & pids+=($!)
done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
"$GRENDEL_EXHAD_PYTHON" "$ALP/merge_manifests.py" merge "$OUT" "$JOBS"
# the manifest over every mass (all bundles exist: nothing is generated; records come from the merged manifest)
# shellcheck disable=SC2086
gen manifest "$OUT" $MASSES || rc=1
n=$(ls "$OUT"/templates_*.npz | wc -l)
echo "$(date '+%F %T') 03_templates done: $n bundles, rc=$rc" >> "$L/03_templates.log"
[ "$rc" -eq 0 ] && [ "$n" -eq 151 ]
