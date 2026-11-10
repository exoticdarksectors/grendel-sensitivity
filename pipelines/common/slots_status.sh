#!/usr/bin/env bash
# Show which slots of the shared CPU pool are taken and by what.
d="${REPRO:-$(cd "$(dirname "$0")/../.." && pwd)}/work/slots"
for f in "$d"/slot*.lock; do [ -s "$f" ] && printf '%s: %s\n' "$(basename "$f" .lock)" "$(cat "$f")"; done
echo "pool: ${GRENDEL_SLOTS:-$(( $(getconf _NPROCESSORS_ONLN) - 1 ))} slots"
q="$d/queue"
if [ -d "$q" ]; then
  n=0; for t in "$q"/[0-9]*; do [ -e "$t" ] || continue; n=$((n+1)); pid="${t##*.}"
    printf 'queue %d: pid %s %s\n' "$n" "$pid" "$(ps -o args= -p "$pid" 2>/dev/null | grep -o -- '--label [^ ]*' | cut -d' ' -f2)"; done
  echo "queued (FIFO): $n"
fi
