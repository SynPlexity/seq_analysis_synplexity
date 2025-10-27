#!/usr/bin/env python3
"""
format_seqs_from_template.py

Reads SynPlexity's Excel order template (tab 2: "Sequences") and outputs:
  - seq_id        (from Excel)
  - seq_to_synth  (cleaned uppercase sequence)

Usage:
    python scripts/format_seqs_from_template.py <input_xlsx>

Behavior:
- Requires columns: 'seq_id' and 'sequence'
- Removes rows where seq_id is empty or only '.' or '…'
- Uppercases and strips whitespace/newlines from sequence
- Writes <basename>.csv to the output/ directory
"""

import sys
import re
from pathlib import Path
import pandas as pd

SEQS_SHEET_NAME = "Sequences"
COL_SEQ_ID = "seq_id"
COL_CORE = "sequence"


def read_sequences_sheet(path_xlsx: str) -> pd.DataFrame:
    try:
        return pd.read_excel(path_xlsx, sheet_name=SEQS_SHEET_NAME)
    except Exception:
        # fallback to second sheet by index if the name isn't present
        return pd.read_excel(path_xlsx, sheet_name=1)


def require_columns(df: pd.DataFrame):
    missing = [c for c in [COL_SEQ_ID, COL_CORE] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")


def clean_seq(s: str) -> str:
    if pd.isna(s):
        return ""
    return re.sub(r"\s+", "", str(s)).upper()


def is_only_dots_or_ellipsis(x) -> bool:
    if pd.isna(x):
        return True
    s = re.sub(r"\s+", "", str(x)).strip()
    s_wo = s.replace(".", "").replace("…", "")
    return s_wo == ""


def make_unique(ids):
    seen = {}
    out = []
    for x in ids:
        base = str(x)
        if base not in seen:
            seen[base] = 1
            out.append(base)
        else:
            seen[base] += 1
            out.append(f"{base}_{seen[base]}")
    return out


def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/format_seqs_from_template.py <input_xlsx>")
        sys.exit(1)

    input_xlsx = Path(sys.argv[1])
    if not input_xlsx.is_file():
        print(f"Error: file '{input_xlsx}' not found.")
        sys.exit(1)

    # Create output directory if needed
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    # Derive output path (output/<basename>.csv)
    output_csv = output_dir / f"{input_xlsx.stem}.csv"

    df = read_sequences_sheet(input_xlsx)
    require_columns(df)

    # Drop placeholder rows (empty / only dots/ellipsis)
    before = len(df)
    df = df[~df[COL_SEQ_ID].apply(is_only_dots_or_ellipsis)].copy()
    dropped = before - len(df)
    if dropped:
        print(f"Dropped {dropped} placeholder row(s) with seq_id of only dots/ellipsis.")

    # Build output DataFrame
    out = pd.DataFrame()
    out["seq_id"] = df[COL_SEQ_ID].astype(str)
    out["seq_to_synth"] = df[COL_CORE].apply(clean_seq)

    # Make seq_id unique if duplicates exist
    if out["seq_id"].duplicated().any():
        out["seq_id"] = make_unique(out["seq_id"])

    # Warn if any empty sequences remain
    n_rows = len(out)
    n_empty = (out["seq_to_synth"] == "").sum()
    if n_empty:
        print(f"Warning: {n_empty}/{n_rows} rows have empty seq_to_synth.")

    out.to_csv(output_csv, index=False)
    print(f"Saved reformatted file to: {output_csv}")
    print(f"Rows written: {n_rows}")


if __name__ == "__main__":
    main()
