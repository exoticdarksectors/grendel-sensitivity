"""Write the FASER2-FPF BC10 curves (SensCalc JSON -> 'mass_GeV  invf_GeV^-1' chains, BNT 1/f = g_Y/246).
python3 -I make_curve_files.py <SensCalc tree> <out dir>"""
import sys, glob, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contour_tools as ct
root, out = sys.argv[1], sys.argv[2]
base = os.path.join(root, 'Sensitivity domains', 'ALP-fermion', 'FASER2-FPF')
for nev, tag in [('3.', 'nev3'), ('2.3', 'nev2p3')]:
    # repro: the glob also matches *.metadata.json; the campaign got the right file from the macOS (APFS) listing order
    src = [f for f in glob.glob(os.path.join(base, 'Sensitivity_ALP-fermion_at_FASER2-FPF_Nev=%s_Npot=*.json' % nev))
           if not f.endswith('.metadata.json')][0]
    chains, modes = ct.load_senscalc_json(src)
    chains = ct.scale(chains, ct.GY_TO_INVF)
    lo, hi = ct.mass_range(chains)
    header = [
        'FASER2 (FPF geometry, SensCalc setup "FASER2-FPF") fermiophilic-ALP (PBC BC10) sensitivity, computed with SensCalc',
        'commit eb2c456 (github.com/maksymovchynnikov/SensCalc, corrected ALP widths: RG sign fix, matched exHad de535a6),',
        'run headless (campaigns/bc10/faser2_senscalc_20261025, see README.md there).',
        'Date: 2026-10-25',
        'Setup FASER2-FPF: decay volume 650-660 m from the IP, 3 m (x) x 1 m (y) box on the beam axis; detector 3 x 1 m,',
        '  26 m long, ECAL; no magnet. SensCalc facility "LHC" (sqrt(s) = 13 TeV spectra; sigma_inel = 72 mb, sigma_bb = 620 mub;',
        '  SensCalc has no 14 TeV option), 3 ab^-1: N_coll = 3000 fb^-1 x 72 mb = 2.16e17.',
        'Selection: >= 2 detectable decay products (opposite-charge pair) reaching the detector; no kinematic cuts.',
        'Settings: decays "Default" (2501.04525 widths, eb2c456 tables); production "Revised (2501.04525)":',
        '  %s;' % ', '.join(modes),
        '  i.e. AP-Central bremsstrahlung and Central Drell-Yan, as the SHiP curve ship_senscalc.dat.',
        'Threshold: N_ev >= %s (zero background).' % nev.rstrip('.') + (' 95% CL convention of FASER2/FPF.' if tag == 'nev3' else ' 90% CL.'),
        'Source y axis g_Y = 2 v/f_a (SensCalc ALP-fermion coupling); stored as 1/f (BNT) = g_Y/246.',
        'Mass range of the region: %.4f-%.4f GeV. Gaps at 0.125-0.14, 0.538-0.555, 0.94-0.974 GeV are SensCalc\'s' % (lo, hi),
        '  pi0/eta/eta\' resonance exclusions. Below ~0.21 GeV the upper edge (1/f = 4.03e-4) is the top of SensCalc\'s',
        '  coupling grid (g_Y = 0.1), not a physical edge (same as in ship_senscalc.dat).',
        'Source file: %s' % os.path.basename(src),
        'Format: chains separated by blank line; each row "mass_GeV  invf_GeV^-1"',
    ]
    path = os.path.join(out, 'faser2_fpf_senscalc_%s.dat' % tag)
    ct.write_dat(path, chains, header)
    print(path, len(chains), 'chains', sum(len(c) for c in chains), 'vertices', 'mass range %.4f-%.4f' % (lo, hi))
