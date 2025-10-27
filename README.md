# SynPlexity Sequence Analysis

A Python-based workflow for **QC**, **feature analysis**, and **diversity measurement** of customer-provided DNA sequences.

> [!NOTE]
> **Before starting:** \
> ***1. Download the `SynPlexity_Order_Template.xlsx` file.*** \
> ***2. Paste your custom sequences into the `Sequences` tab.*** \
> ***3. Save a copy of the completed template in the `seqs/` directory after cloning this repository.*** \
> ***4. Rename the file using a descriptive project name — this name (e.g., `my_project.xlsx`) will be used throughout the workflow in place of `<custom_seqs>.xlsx`.*** \
> ***5. If you’re setting up the workflow for the first time, use the included `test_384_seqs.xlsx` example dataset in the `seqs/` directory to test and verify your setup.***

---

## 🔍 Overview

This pipeline processes and analyzes designed DNA sequences in three reproducible steps:

1. **Input Formatting**  
- Cleans and standardizes SynPlexity’s Excel template (`.xlsx`) into a two-column CSV file.
   - Auto-generates headers if missing
   - Converts all sequences to uppercase DNA format

2. **QC & Analysis**  
- Computes core sequence quality metrics:
   - Detects illegal restriction sites (BspQI, BtsI-v2)
   - Measures length distribution, GC window imbalance, internal repeats, and homopolymer frequency
   - Produces summary `.csv` files and histograms

3. **Similarity Analysis**
 - Quantifies sequence diversity using:
   - **Hamming % identity** (for equal-length sequences)
   - **k-mer Jaccard similarity** (for cross-length comparisons, when applicable)
   - Outputs summary statistics and histograms for visualizing redundancy

---

## 🚀 Installation & Setup

### Clone the Repository

```bash
git clone git@github.com:SynPlexity/seq_analysis_synplexity.git
cd seq_analysis_synplexity
```

Install dependencies if needed:

```bash
pip install pandas numpy matplotlib
```

---

## ▶️ Quick Start

> [!NOTE]
> **Before running the commands below:** \
> - *Replace `<custom_seqs>.xlsx` with your project’s file name.*
> - *Use the included `test_384_seqs.xlsx` example dataset to test and verify the complete workflow.*

### 1. Format Sequences from Template

Place your input Excel file in the `seqs/` directory and run:

```bash
python scripts/1_format_seqs_from_template.py \
    seqs/<custom_seqs>.xlsx

# Example (using test_384_seqs.xlsx)-----------
python scripts/1_format_seqs_from_template.py \
    seqs/test_384_seqs.xlsx
```

- *Creates an `output/` directory (if missing)*
- *Writes a cleaned, standardized CSV to `output/<custom_seqs>.csv`*

### 2. Sequence QC & Feature Analysis

Perform sequence-level QC, restriction-site detection, GC window imbalance analysis, repeat counting, and homopolymer profiling:

```bash
python scripts/2_analyze_seqs.py \
    output/<custom_seqs>.csv

# Example (using test_384_seqs.xlsx)-----------
python scripts/2_analyze_seqs.py \
    output/test_384_seqs.csv
```

- *All QC outputs are written to `output/`*

### 3. Sequence Similarity Analysis

Compute pairwise similarity using both Hamming % identity and k-mer Jaccard:

```bash
python scripts/3_seq_similarity.py \
    -i output/<custom_seqs>.csv

# Example (using test_384_seqs.xlsx)-----------
python scripts/3_seq_similarity.py \
    -i output/test_384_seqs.csv
```

- *Results are written to `output/similarity/`*

---

## 📊 QC & Analysis Outputs

Example output from **Step 2**:

```bash
output/
├── <custom_seqs>_gc_window_outside_35_65.csv   # % of GC windows outside 35–65%
├── <custom_seqs>_longest_repeat.csv            # Longest internal repeat per sequence
├── <custom_seqs>_repeats_per_seq.csv           # Count of internal repeats ≥9 bp per sequence
├── <custom_seqs>_homopolymer_counts.png        # Homopolymer counts per sequence
├── <custom_seqs>_length_distribution.png       # Sequence length histogram
├── <custom_seqs>_gc_window_imbalance.png       # GC window imbalance histogram
├── <custom_seqs>_repeat_length_counts.png      # Total internal repeat frequencies
└── <custom_seqs>.log                           # Run log file
```

---

## 📈 Similarity Analysis Outputs

Example output from **Step 3**:

```bash
output/similarity/
├── summary.txt              # Mean, median, and percentile similarity stats
├── pairwise_sample.csv      # Sampled pairwise similarity values
├── hamming_hist.png         # Histogram of Hamming % identity (equal-length)
├── jaccard_hist.png         # Histogram of k-mer Jaccard similarity (cross-length)
└── combined_hist.png        # Combined similarity distribution
```

**Summary metrics include:**

- **Mean/median similarity:** overall redundancy within the sequence dataset
- **p90 similarity:** top 10% of most-similar pairs (potential duplicates)
- **Histograms:** visualize sequence diversity and clustering of related designs

> [!NOTE]
> The `jaccard_hist.*` files are only generated if the dataset includes sequences of differing lengths.

---

## ⚙️ Maintainers
 
- Karl Romanowicz (karl@synplexity.com)
- Calin Plesa (calin@synplexity.com)