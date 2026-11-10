#!/usr/bin/env bash
# tools/export_public.sh DEST -- write the public version of this repository: the tracked files of HEAD (git archive)
# in a new git repository DEST, as one commit by the configured git user, on branch main. Nothing is pushed.
#
#     bash tools/export_public.sh ../grendel-repro-public
#
# Refuses to run with uncommitted changes, or if a tracked file outside src/ still holds a path of this machine
# (the home directory, /Users/, /Volumes/, /afs/, /eos/). src/ holds verbatim exports of other repositories and is
# left as it is.
set -euo pipefail
REPRO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:?usage: export_public.sh DEST}"
cd "$REPRO"
if [ -e "$DEST" ]; then echo "$DEST exists: choose a new folder" >&2; exit 1; fi
if ! git diff --quiet HEAD -- || [ -n "$(git ls-files --others --exclude-standard)" ]; then
  echo "uncommitted or untracked files: commit them first (git status)" >&2; exit 1
fi
if git grep -nI -e "$HOME" -e '/Users/' -e '/Volumes/' -e '/afs/' -e '/eos/' -- . ':!src' ':!tools/export_public.sh'; then
  echo "the lines above hold paths of this machine: fix them first" >&2; exit 1
fi
mkdir -p "$DEST"
git archive --format=tar HEAD | tar -x -C "$DEST"
cd "$DEST"
git init -q -b main
git add -A
git commit -q -m "Reproduction of the sensitivity figures of the GRENDEL v2 paper (arXiv:2609.00152v2, Sec. IV)"
echo "public tree: $(pwd) ($(git ls-files | wc -l) files, commit $(git rev-parse --short HEAD), no remote)"
