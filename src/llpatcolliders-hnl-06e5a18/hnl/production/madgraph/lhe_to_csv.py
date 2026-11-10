#!/usr/bin/env python3
"""Convert MadGraph LHE files to HNL or tau CSV rows."""

import gzip
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from production.io import write_csv_matrix


class LHEParser:
    PDG_HNL = 9900012
    PDG_TAU_PLUS = -15
    PDG_TAU_MINUS = 15

    def __init__(self, lhe_path):
        self.lhe_path = Path(lhe_path)
        if not self.lhe_path.exists():
            raise FileNotFoundError(f"LHE file not found: {lhe_path}")

    def _open(self):
        if self.lhe_path.suffix == '.gz':
            return gzip.open(self.lhe_path, 'rt', encoding='utf-8')
        return open(self.lhe_path, 'r', encoding='utf-8')

    def extract_particles(self, target_pdg_ids, split_event_weight=False):
        in_event = False
        header_parsed = False
        event_weight = 1.0
        all_event_pdgs = []
        candidate_records = []

        with self._open() as f:
            for line in f:
                stripped = line.strip()

                if stripped.startswith('<event>'):
                    in_event = True
                    header_parsed = False
                    event_weight = 1.0
                    all_event_pdgs = []
                    candidate_records = []
                    continue

                if stripped.startswith('</event>'):
                    if candidate_records:
                        if split_event_weight:
                            particle_weight = event_weight / len(candidate_records)
                        else:
                            particle_weight = event_weight
                        for record, mother_idx in candidate_records:
                            record['weight'] = particle_weight
                            origin = 0
                            if 1 <= mother_idx <= len(all_event_pdgs):
                                origin = all_event_pdgs[mother_idx - 1]
                            record['origin'] = origin
                            yield record
                    in_event = False
                    continue

                if in_event:
                    if stripped.startswith('<') or stripped.startswith('#'):
                        continue

                    if not header_parsed:
                        parts = stripped.split()
                        if len(parts) >= 3:
                            event_weight = float(parts[2])
                        header_parsed = True
                        continue

                    parts = stripped.split()
                    if len(parts) >= 11:
                        try:
                            pdgid = int(parts[0])
                        except (ValueError, IndexError):
                            pdgid = 0
                        all_event_pdgs.append(pdgid)
                        if pdgid in target_pdg_ids:
                            try:
                                mother_idx = int(parts[2])
                                candidate_records.append((
                                    {
                                        'weight': event_weight,
                                        'pdgid': pdgid,
                                        'E': float(parts[9]),
                                        'px': float(parts[6]),
                                        'py': float(parts[7]),
                                        'pz': float(parts[8]),
                                    },
                                    mother_idx,
                                ))
                            except (ValueError, IndexError):
                                pass

    def write_hnl_csv(self, output_path):
        rows = []
        for p in self.extract_particles({self.PDG_HNL}):
            rows.append([p['weight'], p['E'], p['px'], p['py'], p['pz']])
        write_csv_matrix(output_path, rows)
        return len(rows)

    def write_tau_csv(self, output_path):
        rows = []
        for p in self.extract_particles(
            {self.PDG_TAU_MINUS, self.PDG_TAU_PLUS},
            split_event_weight=False,
        ):
            rows.append([p['weight'], p['E'], p['px'], p['py'], p['pz'], p['origin']])
        write_csv_matrix(output_path, rows)
        return len(rows)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Convert LHE to CSV")
    parser.add_argument("lhe_file", help="Input LHE file (.lhe or .lhe.gz)")
    parser.add_argument("csv_file", help="Output CSV file")
    parser.add_argument("--mode", choices=["hnl", "tau"], default="hnl",
                        help="Extraction mode (default: hnl)")
    args = parser.parse_args()

    lhe = LHEParser(args.lhe_file)
    if args.mode == "hnl":
        n = lhe.write_hnl_csv(args.csv_file)
    else:
        n = lhe.write_tau_csv(args.csv_file)

    print(f"Extracted {n} particles → {args.csv_file}")


if __name__ == "__main__":
    main()
