# pipelines/common/ref.sh -- where the validation steps find their reference material. Sourced (bash) by the common
# file of every pipeline, after env/activate.sh.
#
# No production step reads any of it. Three kinds:
#   $GRENDEL_PUBLISHED      the published tables: src/grendel-sensitivity/curves/data (in this repository)
#   $GRENDEL_PAPER_FIGURES  the figure files of the paper: reference/paper_figures (in this repository)
#   $GRENDEL_REF            the outputs of the original runs (campaign records, event samples, templates, ...):
#                           optional and not in this repository (reference/README.md; GRENDEL_REF=none: without
#                           it). Without it every comparison that needs it prints a SKIP line, collected in
#                           work/skipped_validations.txt.
#
#   if have_ref "BC5 four-vectors vs the published pool" "$REF_BC5_VECTORS"; then python validate...; fi
: "${REPRO:?source env/activate.sh first}"
: "${GRENDEL_PUBLISHED:=$REPRO/src/grendel-sensitivity/curves/data}"
: "${GRENDEL_PAPER_FIGURES:=$REPRO/reference/paper_figures}"
REF="${GRENDEL_REF:-none}"
[ "$REF" != none ] || REF=/nonexistent/GRENDEL_REF # so that "$REF/..." never exists when there is no reference tree

# have_ref LABEL PATH... -- true if every PATH exists; otherwise print "SKIP LABEL: ..." and log it, false
have_ref() {
  local label=$1 p n=0 first=""
  shift
  for p in "$@"; do
    if [ ! -e "$p" ]; then n=$((n + 1)); [ -n "$first" ] || first=$p; fi
  done
  [ "$n" -eq 0 ] && return 0
  local what="${first#"$REF"/}"
  [ "$n" -gt 1 ] && what="$what and $((n - 1)) more"
  local msg="SKIP $label: reference data not available ($what; see reference/README.md)"
  echo "$msg" >&2
  mkdir -p "$REPRO/work"
  printf '%s\t%s\n' "$(date '+%F %T')" "$msg" >> "$REPRO/work/skipped_validations.txt"
  return 1
}
