#!/usr/bin/env python3
"""
4_seq_similarity.py

Estimate within-dataset sequence similarity from a SynPlexity-formatted CSV
(two columns: header, seq_to_synth). Uses:
  - Hamming % identity for equal-length pairs (exact, no gaps)
  - k-mer Jaccard similarity for unequal-length pairs (fast proxy)

Outputs under output/similarity/:
  - summary.txt
  - pairwise_sample.csv
  - hamming_hist.png      (if any equal-length pairs)
  - jaccard_hist.png      (if any cross-length pairs)
  - combined_hist.png     (if any values)
  - *_hist_counts.csv     (binned counts for each histogram)

Usage:
  python scripts/4_seq_similarity.py -i seqs/<custom_seqs>.csv
"""

import argparse
import csv
import os
import random
import statistics
from collections import defaultdict
from itertools import combinations

# --- I/O helpers -------------------------------------------------------------

def read_formatted_csv(path):
    seqs = []
    with open(path, newline='') as fh:
        r = csv.DictReader(fh)
        if not r.fieldnames:
            raise ValueError("Input CSV has no header row.")
        cols = {c.lower(): c for c in r.fieldnames}
        header_col = cols.get('header', None)
        seq_col = cols.get('seq_to_synth', cols.get('sequence', None))
        if not seq_col:
            raise ValueError("Could not find 'seq_to_synth' or 'sequence' column in CSV.")
        for row in r:
            header = row[header_col] if header_col else ""
            seq = (row[seq_col] or "").strip().upper()
            if not seq:
                continue
            seqs.append((header, seq))
    if not seqs:
        raise ValueError("No sequences found.")
    return seqs

def ensure_outdir():
    """
    Ensure output/similarity/ exists.
    Returns the absolute path to output/similarity.
    """
    base_outdir = os.path.abspath("output")
    sim_outdir = os.path.join(base_outdir, "similarity")
    os.makedirs(sim_outdir, exist_ok=True)
    return sim_outdir

# --- Similarity metrics ------------------------------------------------------

def hamming_identity(a, b):
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / len(a)

def kmers(s, k):
    if len(s) < k:
        return set()
    return {s[i:i+k] for i in range(len(s)-k+1)}

def jaccard(a_set, b_set):
    if not a_set and not b_set:
        return 1.0
    if not a_set or not b_set:
        return 0.0
    inter = len(a_set & b_set)
    union = len(a_set | b_set)
    return inter / union if union else 0.0

def sample_pairs(n, max_pairs):
    total_pairs = n*(n-1)//2
    if max_pairs >= total_pairs:
        return None  # signal to iterate all pairs
    seen = set()
    while len(seen) < max_pairs:
        i = random.randrange(n)
        j = random.randrange(n)
        if i == j:
            continue
        a, b = (i, j) if i < j else (j, i)
        seen.add((a, b))
    return seen

# --- Summaries & plots -------------------------------------------------------

def summarize(name, values):
    if not values:
        return f"{name}: n=0"
    vals = sorted(values)
    def pct(p):
        idx = max(0, min(len(vals)-1, int(round(p*(len(vals)-1)))))
        return vals[idx]
    return (
        f"{name} (n={len(vals)}):\n"
        f"  mean:   {statistics.mean(vals):.4f}\n"
        f"  median: {statistics.median(vals):.4f}\n"
        f"  p10:    {pct(0.10):.4f}\n"
        f"  p25:    {pct(0.25):.4f}\n"
        f"  p75:    {pct(0.75):.4f}\n"
        f"  p90:    {pct(0.90):.4f}\n"
        f"  min:    {vals[0]:.4f}\n"
        f"  max:    {vals[-1]:.4f}\n"
    )

def plot_hist(values, out_png, out_csv, title, bins=40):
    # Lazy import to keep CLI snappy and avoid importing if not plotting
    import matplotlib
    matplotlib.use("Agg")  # headless
    import matplotlib.pyplot as plt
    import numpy as np

    if not values:
        return False

    # Clamp domain to [0,1] for similarity
    arr = np.clip(np.array(values, dtype=float), 0.0, 1.0)
    counts, edges = np.histogram(arr, bins=bins, range=(0.0, 1.0))

    # Save counts
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bin_left", "bin_right", "count"])
        for i in range(len(counts)):
            w.writerow([f"{edges[i]:.5f}", f"{edges[i+1]:.5f}", int(counts[i])])

    # Plot (single chart, no custom colors/styles)
    plt.figure(figsize=(7, 4))
    plt.hist(arr, bins=bins, range=(0.0, 1.0))
    plt.title(title)
    plt.xlabel("Similarity")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    return True

# --- Main --------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Estimate pairwise sequence similarity.")
    ap.add_argument("-i", "--input", required=True, help="Path to <RUN_NAME>_formatted.csv")
    ap.add_argument("--kmer", type=int, default=5, help="k for k-mer Jaccard (default: 5)")
    ap.add_argument("--sample", type=int, default=50000,
                    help="Max number of total pair comparisons to compute (default: 50000)")
    ap.add_argument("--equal-limit", type=int, default=300,
                    help="If equal-length group size ≤ this, compute all Hamming pairs (default: 300)")
    ap.add_argument("--bins", type=int, default=40, help="Histogram bins (default: 40)")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed for sampling")
    args = ap.parse_args()

    random.seed(args.seed)
    seqs = read_formatted_csv(args.input)

    # <- changed: always write under output/similarity
    outdir = ensure_outdir()

    pairs_csv = os.path.join(outdir, "pairwise_sample.csv")
    summary_txt = os.path.join(outdir, "summary.txt")
    ham_png = os.path.join(outdir, "hamming_hist.png")
    jac_png = os.path.join(outdir, "jaccard_hist.png")
    mix_png = os.path.join(outdir, "combined_hist.png")
    ham_counts = os.path.join(outdir, "hamming_hist_counts.csv")
    jac_counts = os.path.join(outdir, "jaccard_hist_counts.csv")
    mix_counts = os.path.join(outdir, "combined_hist_counts.csv")

    # Group by length and precompute kmers
    by_len = defaultdict(list)
    for hdr, s in seqs:
        by_len[len(s)].append((hdr, s))
    k = args.kmer
    kmaps = [kmers(s, k) for _, s in seqs]

    # Collect similarities
    hamming_vals, jaccard_vals = [], []

    # Write sampled pairs table
    with open(pairs_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id1", "id2", "len1", "len2", "method", "similarity"])

        # Equal-length (Hamming)
        for L, arr in sorted(by_len.items()):
            m = len(arr)
            if m < 2:
                continue
            # do all pairs if small group
            if m <= args.equal_limit:
                idx_iter = combinations(range(m), 2)
            else:
                # sample within this length group
                group_pairs = min(args.sample // max(1, len(by_len)), args.sample)
                idx_pairs = sample_pairs(m, group_pairs)
                idx_iter = combinations(range(m), 2) if idx_pairs is None else idx_pairs
            for i, j in idx_iter:
                id1, s1 = arr[i]
                id2, s2 = arr[j]
                sim = hamming_identity(s1, s2)
                hamming_vals.append(sim)
                w.writerow([id1, id2, L, L, "hamming_identity", f"{sim:.6f}"])

        # Cross-length (Jaccard)
        lengths = [len(s) for _, s in seqs]
        by_len_indices = defaultdict(list)
        for idx, L in enumerate(lengths):
            by_len_indices[L].append(idx)
        lengths_seen = sorted(by_len_indices.keys())
        if len(lengths_seen) > 1:
            remaining = max(0, args.sample - len(hamming_vals))
            for _ in range(remaining):
                L1, L2 = random.sample(lengths_seen, 2)
                i = random.choice(by_len_indices[L1])
                j = random.choice(by_len_indices[L2])
                sim = jaccard(kmaps[i], kmaps[j])
                jaccard_vals.append(sim)
                id1, s1 = seqs[i]
                id2, s2 = seqs[j]
                w.writerow([id1, id2, len(s1), len(s2), "kmer_jaccard", f"{sim:.6f}"])

    # Write summary
    with open(summary_txt, "w") as fh:
        fh.write(f"Input: {args.input}\n")
        fh.write(f"Total sequences: {len(seqs)}\n")
        fh.write(f"k-mer size (Jaccard): {k}\n")
        fh.write(f"Equal-length all-pairs cutoff: {args.equal_limit}\n")
        fh.write(f"Pair sample target: {args.sample}\n\n")
        fh.write("Length distribution (count by length):\n")
        for L in sorted(by_len.keys()):
            fh.write(f"  {L:>5} nt : {len(by_len[L])}\n")
        fh.write("\n")
        fh.write(summarize("Hamming % identity (equal-length)", hamming_vals))
        fh.write("\n")
        fh.write(summarize("k-mer Jaccard (cross-length)", jaccard_vals))
        fh.write("\n")
        all_vals = (hamming_vals or []) + (jaccard_vals or [])
        fh.write(summarize("Combined similarity (mixed metrics)", all_vals))
        fh.write("\nDone.\n")

    # Plots (one chart per file; no styles/colors specified)
    if hamming_vals:
        plot_hist(hamming_vals, ham_png, ham_counts, "Hamming % identity (equal-length)", bins=args.bins)
    if jaccard_vals:
        plot_hist(jaccard_vals, jac_png, jac_counts, "k-mer Jaccard (cross-length)", bins=args.bins)
    if (hamming_vals or jaccard_vals):
        plot_hist(all_vals, mix_png, mix_counts, "Combined similarity", bins=args.bins)

    print(f"[similarity] Wrote: {pairs_csv}")
    print(f"[similarity] Wrote: {summary_txt}")
    if hamming_vals:
        print(f"[similarity] Wrote: {ham_png}, {ham_counts}")
    if jaccard_vals:
        print(f"[similarity] Wrote: {jac_png}, {jac_counts}")
    if (hamming_vals or jaccard_vals):
        print(f"[similarity] Wrote: {mix_png}, {mix_counts}")

if __name__ == "__main__":
    main()
