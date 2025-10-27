#!/usr/bin/env python3
"""
2_analyze_seqs.py

Adds per-sequence internal repeat metrics and logs output, including sequence count statistics:

* repeats_per_seq.csv           — number of perfect internal repeats ≥9 bp per sequence
* longest_repeat_per_seq.csv    — longest repeat length per sequence

Also plots:
    • Histogram of repeats-per-sequence
    • Histogram of longest repeat length per sequence

Prints summary stats for all metrics and the following sequence count bins:
    • Total sequences
    • 0-400 bp
    • 401-600 bp
    • 601-800 bp
    • 801-1000 bp
    • 1001-1200 bp
    • 1201-1400 bp
    • 1401-1600 bp
    • 1601-1800 bp
    • 1801-2000 bp

Previous functionality (restriction-site counts, GC imbalance, homopolymers, etc.)
remains unchanged.

Usage:
    python scripts/2_analyze_seqs.py <custom_seqs>.csv

Dependencies: numpy, pandas, matplotlib
"""

from __future__ import annotations
import sys
import re
from collections import Counter
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ───────────────────────────────────────────────────────────────────────────── #
# Parameters
MIN_REPEAT_LEN = 9
MIN_HOMOPOLYMER = 8
WINDOW = 20
GC_LOW, GC_HIGH = 0.35, 0.65
RE_SITES = {
    "BspQI":  ("GCTCTTC", "GAAGAGC"),
    "BtsI_v2":("GCAGTG",  "CACTGC"),
}
# ───────────────────────────────────────────────────────────────────────────── #

class Tee:
    """Tee stdout and stderr to multiple file-like objects."""
    def __init__(self, *files):
        self.files = files
    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()


def load_dataframe(csv_path: Union[str, Path]) -> pd.DataFrame:
    df = pd.read_csv(csv_path, dtype=str)
    if {"seq_id", "seq_to_synth"} - set(df.columns):
        raise ValueError("CSV must have columns 'seq_id' and 'seq_to_synth'")
    df["seq_to_synth"] = df["seq_to_synth"].str.upper().str.strip()
    return df


def basic_stats(arr: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
    }


def count_internal_repeats(seq: str, min_len: int = MIN_REPEAT_LEN) -> Counter:
    n = len(seq)
    reps = Counter()
    for k in range(min_len, n // 2 + 1):
        seen: dict[str, int] = {}
        for i in range(n - k + 1):
            sub = seq[i:i + k]
            if sub in seen:
                reps[k] += 1
            else:
                seen[sub] = i
    return reps


def homopolymer_stats(seq: str, min_len: int = MIN_HOMOPOLYMER) -> tuple[int, int]:
    pattern = rf"(A{{{min_len},}}|C{{{min_len},}}|G{{{min_len},}}|T{{{min_len},}})"
    lengths = [len(m.group(0)) for m in re.finditer(pattern, seq)]
    return len(lengths), (max(lengths) if lengths else 0)


def pct_windows_gc_outside(seq: str, window: int = WINDOW,
                           low: float = GC_LOW, high: float = GC_HIGH) -> float:
    total = len(seq) - window + 1
    if total <= 0:
        return float('nan')
    low_gc = low * window
    high_gc = high * window
    out = 0
    for i in range(total):
        sub = seq[i:i + window]
        gc = sub.count("G") + sub.count("C")
        if gc < low_gc or gc > high_gc:
            out += 1
    return 100.0 * out / total


def plot_hist(data, bins, xlabel: str, title: str, outfile: Path) -> None:
    plt.figure()
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def main(csv_path: Union[str, Path]) -> None:
    if not csv_path:
        print("Usage: python analyze_seqs_v3.py <input_csv>")
        sys.exit(1)

    input_path = Path(csv_path)
    basename = input_path.stem

    # Create output directory if needed
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    # Write log directly inside output/
    log_path = output_dir / f"{basename}.log"
    log_file = open(log_path, 'w')
    sys.stdout = Tee(sys.stdout, log_file)
    sys.stderr = Tee(sys.stderr, log_file)

    print(f"Logging to {log_path}\n")

    # Load data
    df = load_dataframe(input_path)
    seqs = df["seq_to_synth"].tolist()
    prefix = f"{basename}_"

    # Sequence count stats
    lengths = np.fromiter((len(s) for s in seqs), dtype=int)
    total = len(lengths)
    bins = {
        '0-400': ((lengths <= 400) & (lengths >= 0)).sum(),
        '401-600': ((lengths >= 401) & (lengths <= 600)).sum(),
        '601-800': ((lengths >= 601) & (lengths <= 800)).sum(),
        '801-1000': ((lengths >= 801) & (lengths <= 1000)).sum(),
        '1001-1200': ((lengths >= 1001) & (lengths <= 1200)).sum(),
        '1201-1400': ((lengths >= 1201) & (lengths <= 1400)).sum(),
        '1401-1600': ((lengths >= 1401) & (lengths <= 1600)).sum(),
        '1601-1800': ((lengths >= 1601) & (lengths <= 1800)).sum(),
        '1801-2000': ((lengths >= 1801) & (lengths <= 2000)).sum(),
    }
    print("Sequence count stats:")
    print(f"  Total sequences      : {total}")
    for range_label, count in bins.items():
        print(f"  {range_label:13}: {count}")

    # Restriction site counts
    site_counts = {name: 0 for name in RE_SITES}
    for seq in seqs:
        for name, motifs in RE_SITES.items():
            site_counts[name] += sum(seq.count(m) for m in motifs)
    print("\nRestriction site counts:")
    for name, count in site_counts.items():
        print(f"  {name:8}: {count}")

    # Sequence length stats and plot
    print("\nSequence length stats:", basic_stats(lengths))
    plot_hist(
        lengths, bins=30,
        xlabel="Length (bp)", title="Length distribution",
        outfile=output_dir / f"{prefix}length_distribution.png"
    )

    # GC window imbalance
    pct_out = [pct_windows_gc_outside(s) for s in seqs]
    pct_arr = np.array(pct_out, dtype=float)
    print("\nGC window imbalance % stats:", basic_stats(pct_arr))

    gc_csv = output_dir / f"{prefix}gc_window_outside_35_65.csv"
    pd.DataFrame({"seq_id": df["seq_id"], "pct_gc_out": pct_out}).to_csv(gc_csv, index=False)

    plot_hist(
        pct_arr, bins=30,
        xlabel="% windows outside GC range", title="GC imbalance",
        outfile=output_dir / f"{prefix}gc_window_imbalance.png"
    )

    # Internal repeats
    repeats = []
    longest = []
    agg = Counter()
    for seq in seqs:
        rc = count_internal_repeats(seq)
        repeats.append(sum(rc.values()))
        longest.append(max(rc) if rc else 0)
        agg.update(rc)

    if agg:
        x, y = zip(*sorted(agg.items()))
        plt.figure()
        plt.bar(x, y)
        plt.xlabel("Repeat length (bp)")
        plt.ylabel("Total repeats")
        plt.title("Internal repeats ≥9 bp")
        plt.tight_layout()
        plt.savefig(output_dir / f"{prefix}repeat_length_counts.png", dpi=300)
        plt.close()
    else:
        print("\nNo internal repeats ≥9 bp found.")

    repeats_arr = np.array(repeats)
    longest_arr = np.array(longest)
    print("\nRepeats ≥9 bp per sequence stats:", basic_stats(repeats_arr))
    print("Longest repeat length stats:", basic_stats(longest_arr))

    plot_hist(
        repeats_arr, bins=range(0, repeats_arr.max() + 2),
        xlabel="Repeats per sequence", title="Repeats per sequence",
        outfile=output_dir / f"{prefix}repeats_per_seq_hist.png"
    )
    plot_hist(
        longest_arr, bins=range(0, longest_arr.max() + 2),
        xlabel="Longest repeat (bp)", title="Max repeat per sequence",
        outfile=output_dir / f"{prefix}longest_repeat_hist.png"
    )

    # Save CSVs
    repeats_csv = output_dir / f"{prefix}repeats_per_seq.csv"
    longest_csv = output_dir / f"{prefix}longest_repeat.csv"
    pd.DataFrame({"seq_id": df["seq_id"], "repeats_ge9": repeats}).to_csv(repeats_csv, index=False)
    pd.DataFrame({"seq_id": df["seq_id"], "longest_repeat": longest}).to_csv(longest_csv, index=False)

    # Homopolymers
    hp_counts, hp_max = [], []
    for seq in seqs:
        c, m = homopolymer_stats(seq)
        hp_counts.append(c)
        hp_max.append(m)
    print("\nHomopolymer length stats:", basic_stats(np.array(hp_max)))

    plot_hist(
        hp_counts, bins=range(0, max(hp_counts) + 2),
        xlabel="Homopolymers per sequence", title="Homopolymer counts",
        outfile=output_dir / f"{prefix}homopolymer_counts.png"
    )
    plot_hist(
        np.array(hp_max),
        bins=range(MIN_HOMOPOLYMER, max(hp_max) + 2),
        xlabel="Max homopolymer (bp)", title="Max homopolymer length",
        outfile=output_dir / f"{prefix}homopolymer_max.png"
    )

    print(f"\nAnalyses complete. Outputs in {output_dir}, log in {log_path}.")
    log_file.close()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) == 2 else None)
