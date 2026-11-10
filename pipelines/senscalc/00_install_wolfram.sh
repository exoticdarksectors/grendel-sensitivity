#!/usr/bin/env bash
# 00_install_wolfram.sh -- install the free Wolfram Engine for developers (Linux x86-64) in user space.
#
#   source env/activate.sh && bash pipelines/senscalc/00_install_wolfram.sh
#
# Downloads the official installer (Wolfram Engine 15.0.0, the version the reference campaigns used on macOS)
# to work/senscalc/downloads/ (resumable) and installs it non-interactively with
#     WolframEngine_15_LIN.sh -- -auto -wolframscript=n -targetdir=$REPRO/env/wolfram/WolframEngine -execdir=$REPRO/env/wolfram/bin
# so that env/activate.sh puts env/wolfram/bin (wolframscript, WolframKernel, wolfram, math) on PATH.
# The installer refuses non-root installs unless -auto -wolframscript=n is given (MathInstaller, InstallAsRoot()):
# the WolframScript .deb it would install system-wide is unpacked instead with `dpkg-deb -x` into
# env/wolfram/WolframScript, and env/wolfram/bin/wolframscript is a two-line wrapper that points it at our kernel
# (WOLFRAMSCRIPT_KERNELPATH, unless already set).
# No sudo, nothing outside $REPRO except the licence that `wolframscript -activate` writes later to
# ~/.WolframEngine/Licensing (the activation itself needs the user's Wolfram ID and is NOT done here).
#
# Resume-safe: the download continues where it stopped; the install is skipped when the kernel already exists.
set -euo pipefail
: "${REPRO:?source env/activate.sh first}"

VERSION=15.0.0
URL="https://account.wolfram.com/dl/WolframEngine?version=${VERSION}&platform=Linux"
DL="$REPRO/work/senscalc/downloads"
INST="$DL/WolframEngine_15_LIN.sh"
SIZE=2306691927                       # content-length of the 15.0.0.0 installer on files.wolframcdn.com (2026-06-02)
TARGET="$REPRO/env/wolfram/WolframEngine"
EXECDIR="$REPRO/env/wolfram/bin"
mkdir -p "$DL" "$REPRO/env/wolfram"

# ---- download (resumable; the account.wolfram.com link redirects to a signed CDN URL) ----
have=$( [ -f "$INST" ] && stat -c %s "$INST" || echo 0 )
if [ "$have" -ne "$SIZE" ]; then
  echo "downloading Wolfram Engine $VERSION installer ($have / $SIZE bytes present)"
  for attempt in 1 2 3 4 5; do
    # signed CDN links expire: re-resolve the redirect on every attempt, then resume with -C -
    cdn=$(curl -sI --max-time 60 "$URL" | awk 'tolower($1)=="location:"{print $2}' | tr -d '\r')
    [ -n "$cdn" ] || { echo "could not resolve the download link"; sleep 30; continue; }
    curl -fL --retry 5 --retry-delay 10 -C - -o "$INST" "$cdn" && break || { echo "attempt $attempt failed"; sleep 30; }
  done
  have=$(stat -c %s "$INST")
  [ "$have" -eq "$SIZE" ] || { echo "download incomplete: $have / $SIZE"; exit 1; }
fi
[ -s "$DL/WolframEngine_15_LIN.sh.sha256" ] || sha256sum "$INST" > "$DL/WolframEngine_15_LIN.sh.sha256"; cat "$DL/WolframEngine_15_LIN.sh.sha256"

# ---- install ----
KERNEL="$TARGET/Executables/WolframKernel"
if [ -x "$KERNEL" ]; then
  echo "Wolfram Engine already installed in $TARGET"
else
  rm -rf "$TARGET"; mkdir -p "$EXECDIR"
  # -auto: no questions; -wolframscript=n: no system package (allows a non-root install);
  # -targetdir/-execdir: user-space layout; -nodesktop: no menu entries / icons
  bash "$INST" -- -auto -wolframscript=n -targetdir="$TARGET" -execdir="$EXECDIR" -nodesktop -verbose
fi
[ -x "$KERNEL" ] || { echo "ERROR: no $KERNEL after the install"; exit 1; }

# ---- wolframscript: unpack the bundled .deb in user space and wrap it ----
WSDIR="$REPRO/env/wolfram/WolframScript"
if [ ! -x "$WSDIR/usr/bin/wolframscript" ] && [ ! -x "$WSDIR/opt/Wolfram/WolframScript/bin/wolframscript" ]; then
  deb=$(ls "$TARGET"/SystemFiles/Installation/*_amd64.deb | head -1)
  echo "unpacking $deb"
  rm -rf "$WSDIR"; mkdir -p "$WSDIR"
  dpkg-deb -x "$deb" "$WSDIR"
fi
WSBIN=$(find "$WSDIR" -type f -name wolframscript -perm -u+x | head -1)
[ -n "$WSBIN" ] || { echo "ERROR: no wolframscript binary in $WSDIR"; exit 1; }
rm -f "$EXECDIR/wolframscript"
cat > "$EXECDIR/wolframscript" <<WRAP
#!/bin/sh
# wrapper written by pipelines/senscalc/00_install_wolfram.sh: user-space WolframScript bound to env/wolfram's kernel
: "\${WOLFRAMSCRIPT_KERNELPATH:=$KERNEL}"
export WOLFRAMSCRIPT_KERNELPATH
exec "$WSBIN" "\$@"
WRAP
chmod +x "$EXECDIR/wolframscript"
ls -l "$EXECDIR"
"$EXECDIR/wolframscript" -version || true
du -sh "$TARGET"
echo "installed. Next: wolframscript -activate (needs the user's Wolfram ID), then 02_check_engine.sh"
