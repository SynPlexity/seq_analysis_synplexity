#!/usr/bin/env python3
"""
csv_to_fasta.py

Convert a two-column CSV of sequences into FASTA format and save to a file in the "output/fasta" directory.

Usage:
    python csv_to_fasta.py <input_csv>

Assumes the input CSV has columns:
    seq_id, seq_to_synth
Each row is written as:
    >seq_id
    seq_to_synth
"""

import sys
from pathlib import Path
import pandas as pd


def main(input_csv: str) -> None:
    # Validate argument
    if not input_csv:
        print("Usage: python csv_to_fasta.py <input_csv>", file=sys.stderr)
        sys.exit(1)

    in_path = Path(input_csv)
    if not in_path.is_file():
        print(f"Error: '{in_path}' not found or is not a file.", file=sys.stderr)
        sys.exit(1)

    # Prepare output directory
    fasta_dir = Path('output') / 'fasta'
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # Read CSV
    df = pd.read_csv(in_path, dtype=str)
    required = {'seq_id', 'seq_to_synth'}
    if not required.issubset(df.columns):
        print(f"Error: Input CSV must contain columns: {required}", file=sys.stderr)
        sys.exit(1)

    # Derive output path
    base = in_path.stem
    out_path = fasta_dir / f"{base}.fasta"

    # Write FASTA
    with open(out_path, 'w') as fasta_file:
        for _, row in df.iterrows():
            header = row['seq_id']
            seq = row['seq_to_synth']
            fasta_file.write(f">{header}\n")
            fasta_file.write(f"{seq}\n")

    print(f"FASTA file written to: {out_path}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
