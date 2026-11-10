#!/usr/bin/env bash
# 02_check_engine.sh -- after `wolframscript -activate`: version, licence limits, C compiler, FeynCalc.
#   source env/activate.sh && bash pipelines/senscalc/02_check_engine.sh
# Output: logs/senscalc/02_check_engine.log. Fails (exit 3) if the engine is not activated.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_kernel
with_slots 1 02_check \
  wolframscript -file "$WL/check_soft.wls" < /dev/null 2>&1 | tee "$SCL/02_check_engine.log"
grep -q 'Compile C test: 4.' "$SCL/02_check_engine.log" || { echo "ERROR: C compilation does not work"; exit 1; }
grep -q 'version: 10.2.1' "$SCL/02_check_engine.log" || { echo "ERROR: FeynCalc 10.2.1 does not load"; exit 1; }
echo "engine OK"
