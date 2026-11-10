"""Make a SensCalc tree behave on a case-sensitive file system (Linux ext4) as on macOS APFS (case-insensitive).

    python3 -I case_links.py <SensCalc tree>            # report case-mismatched path literals
    python3 -I case_links.py <SensCalc tree> --apply    # also create the aliases below

Root cause of the first Linux run (FASER2-FPF upper edge -17 %, SHiP -13 %): two SensCalc notebooks build paths with
a lower-case directory name,
    codes/LLP distribution/prod-pheno-bremsstrahlung-proton.nb:
        gppA[mLLP_, model] = Import[FileNameJoin[{..., "phenomenology", model, "production probabilities",
                                    "Coefficients-" <> model <> "-Lambda=1.-TeV.mx"}], "MX"][[-1, 2]]
    codes/LLP distribution/prod-pheno-ALP-gluon.nb: ".../production probabilities/BrRatios-Msquared-LightMesonDecays-..."
while the directory is "phenomenology/<model>/Production probabilities". On macOS the path resolves; on Linux Import
returns $Failed silently, gppA (the ALP-nucleon coupling, with the pi0/eta/eta' mixing poles) drops out of
FormFactorSquaredBrem, and the bremsstrahlung production probability notebook 2 writes
(phenomenology/ALP-fermion/Production probabilities/ProdProb_ALP-fermion_Bremsstrahlung-AP_<facility>.m) is too
small by |gppA|^2 (x 29 at 0.2 GeV, x 704 at 1 GeV).

Fix, without touching SensCalc code: in every phenomenology/<model>/ that has "Production probabilities", a relative
symlink "production probabilities" -> "Production probabilities". No SensCalc code lists the contents of
phenomenology/<model>/ (checked: FileNames is only used on other directories), so the alias adds nothing else.
The report mode lists every string literal of the .nb/.wl files that matches an on-disk name only up to case;
the others found in eb2c456 are cell/section titles or 4. Plots.nb outputs, not paths read by the BC10 chain.
"""
import os
import re
import sys

LIT = re.compile(r'\\"\\<(.*?)\\>\\"|"((?:[^"\\]|\\.)*)"')
ALIASES = {"production probabilities": "Production probabilities"}  # literal -> name on disk, in phenomenology/<model>/


def literals(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    out = set()
    for m in LIT.finditer(s):
        t = m.group(1) if m.group(1) is not None else m.group(2)
        t = t.replace('\\<', '').replace('\\>', '').strip('"\\ ')
        for comp in t.split("/"):
            comp = comp.strip()
            if 2 < len(comp) < 120:
                out.add(comp)
    return out


def walk(root):
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x != ".git"]
        yield d, dirs, files


def main() -> int:
    root, apply = sys.argv[1], "--apply" in sys.argv
    lits, names = set(), {}
    for d, dirs, files in walk(root):
        for f in files:
            if f.endswith((".nb", ".wl", ".wls")) and os.path.getsize(os.path.join(d, f)) < 20_000_000:
                lits |= literals(os.path.join(d, f))
        for n in dirs + files:
            names.setdefault(n.lower(), set()).add(n)
    hits = [(lit, sorted(names[lit.lower()])) for lit in sorted(lits)
            if lit.lower() in names and lit not in names[lit.lower()]]
    for lit, cands in hits:
        print(f"case-mismatched literal {lit!r}: on disk {cands}{'  -> aliased' if lit in ALIASES else ''}")
    made = 0
    if apply:
        ph = os.path.join(root, "phenomenology")
        for model in sorted(os.listdir(ph)):
            for alias, actual in ALIASES.items():
                d = os.path.join(ph, model)
                if os.path.isdir(os.path.join(d, actual)) and not os.path.lexists(os.path.join(d, alias)):
                    os.symlink(actual, os.path.join(d, alias))
                    made += 1
    print(f"{len(hits)} case-mismatched literals; {made} aliases created")
    return 0


if __name__ == "__main__":
    sys.exit(main())
