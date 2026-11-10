"""Compare the deterministic part of inspect_acceptance.wls output (ours) with the campaign's logs/inspect_acceptances.txt.

    python3 -I compare_inspect.py <ours.txt> <campaign inspect_acceptances.txt> <experiment>

Compared line by line (the "file:" path line excluded): facility/ECAL/table dimensions, particle numbers, decay
description, acceptance mass grid, Br_vis(m) and ctau(m, y=1). These come from the decay tables and the facility
definitions, not from the Monte Carlo, so they must agree exactly (numbers to 1e-12 relative).
"""
import re
import sys


def block(text, exp):
    out, on = [], False
    for line in text.splitlines():
        if line.startswith("file: "):
            on = f"/Acceptances/{exp}/" in line
            continue
        if on:
            out.append(line.rstrip())
    return out


NUM = re.compile(r"[-+]?\d+\.?\d*(?:\*\^[-+]?\d+)?")


def nums(line):
    return [float(x.replace("*^", "e")) for x in NUM.findall(line)]


def main() -> int:
    ours, ref, exp = open(sys.argv[1]).read(), open(sys.argv[2]).read(), sys.argv[3]
    a, b = block(ours, exp), block(ref, exp)
    if not a:
        a = [l.rstrip() for l in ours.splitlines() if not l.startswith("file: ") and not l.startswith("SCEnv|")]
    bad = 0
    for la, lb in zip(a, b):
        if la == lb:
            continue
        na, nb = nums(la), nums(lb)
        if len(na) == len(nb) and re.sub(NUM, "#", la) == re.sub(NUM, "#", lb) and all(
                abs(x - y) <= 1e-12 * max(abs(x), abs(y), 1e-300) for x, y in zip(na, nb)):
            continue
        bad += 1
        print(f"  DIFF ours: {la[:200]}\n       ref : {lb[:200]}")
    if len(a) != len(b):
        bad += 1
        print(f"  line count ours {len(a)} ref {len(b)}")
    print(f"{exp}: acceptance metadata {'IDENTICAL' if bad == 0 else f'{bad} differing lines'} ({len(b)} lines)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
