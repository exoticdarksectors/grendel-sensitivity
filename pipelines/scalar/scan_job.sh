#!/usr/bin/env bash
# Runs one scan command inside its CPU slot (called by scan_one.sh):   scan_job.sh OUT_DIR CMD...
# Memory gate: the jobs start one at a time through $W/.memgate, each when >= 3 GB are available, and the gate is
# held for 90 s after a start so the next job sees this one's memory (a BC5 job reaches ~2.5 GB in its first minute).
# Writes OUT_DIR/scan.log (with /usr/bin/time -v) and, on success, OUT_DIR/time.txt.
O=$1; shift
W="$REPRO/work/scalar"
exec 8>"$W/.memgate"
flock 8
while [ "$(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo)" -lt 3 ]; do sleep 20; done
t0=$(date +%s)
/usr/bin/time -v "$@" > "$O/scan.log" 2>&1 &
pid=$!
for _ in $(seq 18); do kill -0 $pid 2>/dev/null || break; sleep 5; done
flock -u 8; exec 8>&-
wait $pid || { echo "FAIL $(tail -n 3 "$O/scan.log")"; exit 3; }
echo "$(( $(date +%s) - t0 )) s wall, max RSS $(awk '/Maximum resident/ {print int($6/1024)}' "$O/scan.log") MB, host $(hostname)" > "$O/time.txt"
