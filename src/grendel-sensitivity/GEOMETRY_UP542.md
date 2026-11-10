# The GRENDEL tunnel (UP542) from the GRENDEL geometry study — geometry used since October 2026

**CERN-internal** (derived from the CERN civil drawings of the gallery and the restricted 2015 laser scan by the CMS
survey team): this repository is private for that reason. The geometry is kept in **one place**, the GRENDEL geometry study
(https://gitlab.cern.ch/grendel/grendel-geometry, `results/geometry/`, CERN login); this repository mirrors it
unchanged in `grendel/geometry/data/up542/` (checked against its `manifest.json` on import; update with
`python scripts/sync_up542_geometry.py <clone of the study>`; check that it is the latest GitLab version with
`python scripts/sync_up542_geometry.py --check <clone of the study>`, which needs a CERN login and exits 1 if not). grendelSim and llpatcolliders mirror the same files.
Never edit the mirror: change the study and re-sync.

Mirrored release: **grendel-geometry v1.0** (git tag `v1.0`, 16 October 2026; geometry version 2026-10-16.1), the
final release of the study.

## Names used here

| Name | What it is |
|---|---|
| CERN civil drawings (2004) | the civil-engineering drawings of the gallery (structure 3529, LHC-DGS135294001/2/3/10), drawn by the contractor (Dragados–SELI joint venture) and re-issued after construction with the status “as built” on 20 Dec 2004: alignment, floor levels and gradients as built; the cross-section kept its design dimensions (2.60 × 2.774 m) |
| 2015 laser scan | the point cloud of the gallery (CH 56–109) by the CMS survey team, report CMS-I-UR-0261 v3 (EDMS 1568271), data EDMS 1568275; re-registered in 2017; the rotation that 2017 fit left free corrected by the geometry study |
| CERN 3D civil model | CERN's CAD model of the Point 5 caverns, shafts and galleries (PLM ST0116966_01), a design model; gives the rock and concrete |
| grendel-geometry | the GRENDEL geometry study (gitlab.cern.ch/grendel/grendel-geometry): centreline and floor from the 2004 drawings, section from the 2015 scan; the files mirrored here |
| arXiv:2609.00152 tunnel | the tunnel of releases v1–v2.2 of this repository, used for the paper's curves |

## The tunnel of arXiv:2609.00152 (releases v1–v2.2)

Releases v1–v2.2 (the curves of arXiv:2609.00152v1/v2) used a 47-point centreline at a fixed Y = 22 m with a
2.90 m × 3.15 m half-ellipse section in a level frame. It is replaced by the grendel-geometry tunnel; the old one is in
the git tags `v1` … `v2.2`.

## What changed with respect to release v2.2

| | v2.2 (arXiv:2609.00152) | grendel-geometry v1.0 (geometry version 2026-10-16.1) |
|---|---|---|
| Height | tunnel centroid at Y = 22 m everywhere | finished floor +23.80 m (PX56 end) → +22.19 m (PM54 junction, CH 73.06) → +22.81 m (south end), centroid 1.345 m above it (2004 drawings; the 2015 scan agrees to a few cm once the one rotation its 2017 re-fit left free is fixed): the old tunnel was **1.5–3.1 m too low** |
| Plan | 47 points, a **mirror image** of the gallery, 1.7 m short of the PX56 face and 1.3 m into the rock at the end | alignment of the 2004 drawings (straight legs and arcs), CH −1.255 (PX56 shaft face) to 108.505 (end wall), **109.76 m**; z<sub>CMS</sub> = −X<sub>b</sub> (the PX56 end at z ≈ −14 m) |
| Section | 2.90 × 3.15 m, half-ellipse roof | **from the 2015 laser scan** (CMS survey team; 66 sections, CH 55–109): floor 2.88 m, **2.98 m wide** at 0.8 m, springline 1.654 m, crown **2.946 m**, circular roof R 1.517 m; area 7.80 m² (2004 drawings, as designed: 2.60 × 2.774 m, R 1.30) |
| Frame | level | the **official CMS physicist frame**: the level survey frame pitched by **12.34 mrad = 0.707°** (CMS-I-UR-0281); section "up" = the true vertical |
| From IP5 | 24.2–54.9 m, \|η\| up to 0.52 | section centres 25.8–54.6 m, \|η\| ≤ 0.50 |

The curves of tag `v2.3` (BC4, the HNLs, BC5; BC10 with the v2.2 ALP widths, which are being revised) use the
grendel-geometry tunnel: it is higher and slightly farther, so the acceptance drops by about 10 %.

## Frames

- Level survey frame (the frame of the CERN 3D civil model = the CMS survey frame to ~1 cm): origin IP5, X<sub>b</sub> along the
  beam (+ towards PX56 / R56), Y<sub>t</sub> horizontal (+ away from USC55), Z up.
- CMS physicist frame (this code): x = −Y<sub>t</sub>, y = Z cos a + X<sub>b</sub> sin a, z = −X<sub>b</sub> cos a + Z sin a,
  a = 12.34 mrad; the beam rises towards L54.
- The 2017 "Génie Civil" coordinates of the 2015 scan and of the 2017 survey marks in the gallery are CMS physicist coordinates
  with relabelled axes: X<sub>gc</sub> = X<sub>ph</sub>, Y<sub>gc</sub> = −Z<sub>ph</sub> (along the beamline, + towards R56),
  Z<sub>gc</sub> = Y<sub>ph</sub>, as defined in the scan report CMS-I-UR-0261 v3 (EDMS 1568271), corrected by the CMS survey team on
  13 October 2026.

## Material between IP5 and the gallery (for background studies)

`grendel/geometry/data/up542/UP542_material_budget.csv`: rock and concrete (m) and interaction lengths on the straight
line from IP5 to each section centre, every 0.5 m; constants in `UP542_material.json`.

- **21–116 λ<sub>I</sub>** along the gallery, with two thin places: **21.4 ± 0.6** at the PX56 mouth (CH −1, where the line climbs the
  open shaft; 56 by CH 4) and **39.9 ± 0.5** at CH 56, above CMS's cryogenic chimney; up to 115.7 ± 1.8 λ<sub>I</sub> at the south end
  (estimated 1σ: molasse and concrete density, CMS, magnet and gallery geometry in quadrature, 1.1–2.7 % of the total;
  Figure U10 of the study's report). No fiducial cut is applied: the fiducial volume is still to be
  decided. Away from the PX56 mouth the thinnest rock is 5.9 m (rock + concrete 7.7 m), at CH 56. Since version
  2026-10-14.3 the rock and concrete are computed for the section from the 2015 scan on the floor of the 2004 drawings
  (the gallery of the CERN 3D civil model was up to 16 cm too high on the west leg).
- **Molasse 2.50 ± 0.05 g/cm³**, wet bulk density of Point 5 core samples (GADZ report 7223-2, EDMS 1952965), λ<sub>I</sub> from its
  composition with pore water: **λ = 38.3 cm** (standard rock in older estimates: 38.2–41.5 cm). Concrete 2.4 ± 0.1 g/cm³,
  assumed (UXC55 vault lining 1.1–1.5 m on the line, the pillar up to 9 m on the south leg, the gallery lining); the 2004
  drawings give strength classes only (C25/30–C28/35, ordinary concrete), no density.
- **CMS 15–24 λ<sub>I</sub>** (median 20; single lines dip where they cross a crack between HCAL barrel wedges or a gap in the muon system), from the CMS simulation geometry (CMSSW_20_1_0_pre2, the latest tag of master) with the solenoid replaced by the real magnet:
  CMSSW models the cold mass as about 445 t of a copper/liquid-helium mix (12.8 X<sub>0</sub>) and cuts the cryogenic chimney
  through the coil; the real cold mass is 220 t of aluminium-stabilised conductor, 3.9 X<sub>0</sub> (JINST 3 S08004, Table 2.1;
  222.5 t as built, 89.5 % aluminium: MTB93 magnet-test documentation, EDMS 753430; radii: Kircher et al., IEEE Trans. Appl.
  Supercond. 10 (2000) 407), and the cryogenic chimney is a 450 mm vacuum vessel welded to the outer vacuum tank, so the coil is
  not interrupted (Levesy et al., IEEE Trans. Appl. Supercond. 16 (2006) 517). Both points confirmed by the CMS magnet group
  (15 Oct 2026). Magnet 1.58 λ<sub>I</sub> on most lines (CMSSW 2.09), 1.41–1.61 through the opening (CMSSW 0), the chimney's
  walls, shield and pipes included (from the CMS magnet group's chimney drawings; details in the study's report).
- λ<sub>I</sub> is the nucleon interaction length; K<sub>L</sub> interaction lengths are longer.
- milliQan (slab detector at CH 88–93) is not included: it may be removed.
- The PM54 branch (`UP542_pm54_branch.csv`, axis from the walls in the 2015 scan): 56–65 λ<sub>I</sub>
  to its section centres.

## Sources

Public: CMS-I-UR-0261 v3 (EDMS 1568271, the 2015 scan report), CMS-I-UR-0278 (EDMS 1817542), CMS-I-UR-0281
(EDMS 2028417, the CMS frame), CMS-I-UR-0299 (EDMS 2875513, milliQan slab support), CMS JINST 3 S08004, CMS Magnet TDR (CERN/LHCC 97-10), Levesy et al. (IEEE Trans. Appl. Supercond. 16 (2006) 517, EDMS 582236), Kircher et al. (IEEE Trans. Appl. Supercond. 10 (2000) 407), MTB93 minutes (EDMS 753430), CMSSW_20_1_0_pre2, the latest tag of master,
Fern et al., Tunn. Undergr. Space Technol. 77 (2018) 249, PDG. CERN-internal: CERN civil drawings LHC-DGS135294001/2/3
(structure 3529, “as built” 2004), scan data EDMS 1568275, CERN 3D civil model PLM ST0116966_01, HL-LHC Point 5 geotechnical documents
(EDMS 1952965, 1684245, 1763507), the CMS cryogenic chimney drawings CMS-SC-ED-0100/-0107 (EDMS 803355/803380). Full account: the study's report (`report/REPORT.pdf` in the GitLab repository).
