# UP542 geometry for the GRENDEL code repositories (grendel-geometry release v1.0, geometry version 2026-10-16.1)

The single source of the UP542 (GRENDEL tunnel) geometry. Kept in the CERN GitLab repository
https://gitlab.cern.ch/grendel/grendel-geometry (`results/geometry/`, made by `code/scripts/export_geometry.py`);
grendel-sensitivity, grendelSim and llpatcolliders mirror this folder unchanged and check it against `manifest.json`.
**CERN-internal** (derived from the CERN civil drawings of UP542 and the restricted 2015 laser scan by the CMS survey team):
mirror only into private repositories. Change the geometry in the study, never in a mirror.

Sources:
- **CERN civil drawings (2004)**: the drawings of structure 3529 (LHC-DGS135294001/2/3/10), which the contractor
  (Dragados–SELI joint venture) re-issued after construction with the status “as built” on 20 Dec 2004: alignment, floor
  levels and gradients as built; the cross-section kept its design dimensions. Plan and floor here come from our
  digitisation of them (10 Oct 2026).
- **2015 laser scan (CMS survey team)**: CH 56–109, report CMS-I-UR-0261 v3 (EDMS 1568271), data EDMS 1568275, re-registered
  in 2017; the rotation the 2017 re-fit left free (13.0 mrad) fixed by this study. It gives the cross-section used here.
- **CERN 3D civil model** (PLM ST0116966_01): rock and concrete on the lines of sight. **CMS simulation geometry**
  (CMSSW_20_1_0_pre2) with the real magnet: CMS material.
- This geometry replaces the tunnel of the GRENDEL paper, arXiv:2609.00152 (grendel-sensitivity v2.2): 2.90 × 3.15 m,
  level, 1.5–3.1 m too low.

| File | Content |
|---|---|
| `UP542_centreline.csv` | every ~0.24 m of chainage: plan (Xb, Yt) and finished floor level from the 2004 drawings, section centroid, and the centroid in the CMS physicist frame (x, y, z) |
| `UP542_alignment_elements.csv` | the analytic alignment of the 2004 drawings: straight legs and circular arcs (CH −1.255 to 108.505, 109.76 m) |
| `UP542_section.json`, `UP542_section_profile.csv` | section from the 2015 laser scan (JSON key `measured`): floor 2.884 m, 2.975 m wide at 0.8 m, springline 1.654 m, crown 2.946 m, circular roof R 1.517 m, area 7.80 m², centroid 1.345 m above the floor; the section of the 2004 drawings (key `design`: 2.60 × 2.774 m, R 1.30) for reference |
| `UP542_frame.json` | frames: the level survey frame, the CMS physicist frame (pitch 12.34 mrad = 0.707°, CMS-I-UR-0281), the scan's Génie Civil labels (CMS-I-UR-0261 v3) |
| `UP542_landmarks.json` | ends, the red door (CH 50.9), the PM54 junction (CH 73.06), floor levels printed on the 2004 drawings, milliQan |
| `UP542_pm54_branch.csv` | the branch from the PM54 junction to the PM54 shaft: axis from the walls in the 2015 scan, floor from the levels of the 2004 drawings; section of the 2004 drawings (2.60 m wide) |
| `UP542_material_budget.csv`, `UP542_material.json` | material from IP5 to each section centre (metres and λ<sub>I</sub>): molasse 2.50 ± 0.05 g/cm³ (Point 5 core samples, GADZ report 7223-2), λ 38.3 cm; concrete 2.40 ± 0.10 g/cm³ (assumed); CMS from the CMS simulation geometry (CMSSW_20_1_0_pre2) with the real magnet |

Conventions: the gallery centroid follows the floor (+23.80 m at PX56, +22.19 m at the PM54 junction, +22.81 m at the
south end, plus 1.345 m); sections are perpendicular to the centreline with "up" = the true vertical, which is
(0, cos a, sin a) in CMS physicist coordinates. Handedness: z_CMS = −Xb (the PX56 end is at z_CMS ≈ −14 m).
