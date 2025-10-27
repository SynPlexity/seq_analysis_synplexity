# SynPlexity Sequence Analysis

A Python-based workflow for **QC**, **feature analysis**, and **diversity measurement** of customer-provided DNA sequence libraries.

> [!NOTE]
> **Before starting:**
>   **1. Download the SynPlexity_Order_Template.xlsx file.**
>   **2. Paste your custom sequences into the Sequences tab.**
>   **3. Save a copy of the completed template in the seqs/ directory after cloning this repository.**
>   **4. Rename the file using a descriptive project name — this name (e.g., `my_project.xlsx`) will be used throughout the workflow in place of `<custom_seqs>.xlsx`.**

---

## 🔍 Overview

This pipeline processes and analyzes designed DNA libraries in four reproducible steps:

1. **Input Formatting**  
- Cleans and standardizes SynPlexity’s Excel template (`.xlsx`) into a two-column CSV file.
   - Auto-generates headers if missing
   - Converts all sequences to uppercase DNA format

2. **FASTA Export**  
- Converts the formatted CSV into FASTA format for reference and downstream tools.

3. **QC & Analysis**  
- Computes core sequence quality metrics:
   - Detects illegal restriction sites (BspQI, BtsI-v2)
   - Measures length distribution, GC window imbalance, internal repeats, and homopolymer frequency
   - Produces summary `.csv` files and histograms

4. **Similarity Analysis**
 - Quantifies library diversity using:
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

### 1. Format Sequences from Template

Place your input Excel file in the `seqs/` directory and run:

> [!NOTE]
> Replace `<custom_seqs>` with a relevant file name for your project

```bash
python scripts/1_format_seqs_from_template.py \
    seqs/<custom_seqs>.xlsx
```

- Creates an `output/` directory (if missing)
- Writes a cleaned, standardized CSV:
    ```bash
    output/<custom_seqs>.csv
    ```

### 2. Export CSV as FASTA

Convert the formatted CSV into FASTA format for downstream use:

```bash
python scripts/2_csv_to_fasta.py \
    output/<custom_seqs>.csv
```

- Creates a subdirectory:
    ```bash
    output/fasta/<custom_seqs>.fasta
    ```

### 3. Sequence QC & Feature Analysis

Perform sequence-level QC, restriction-site detection, GC window imbalance analysis, repeat counting, and homopolymer profiling:

```bash
python scripts/3_analyze_seqs.py \
    output/<custom_seqs>.csv
```

- All QC outputs are written to:
    ```bash
    output/
    ```

### 4. Sequence Similarity Analysis

Compute pairwise similarity using both Hamming % identity and k-mer Jaccard:

```bash
python scripts/4_seq_similarity.py \
    -i output/<custom_seqs>.csv
```

- Results are written to:
    ```bash
    output/similarity/
    ```

---

## 📊 QC & Analysis Outputs

Example output from **Step 3**:

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

Example output from **Step 4**:

```bash
output/similarity/
├── summary.txt              # Mean, median, and percentile similarity stats
├── pairwise_sample.csv      # Sampled pairwise similarity values
├── hamming_hist.png         # Histogram of Hamming % identity (equal-length)
├── jaccard_hist.png         # Histogram of k-mer Jaccard similarity (cross-length)
└── combined_hist.png        # Combined similarity distribution
```

**Summary metrics include:**

- **Mean/median similarity:** overall redundancy within the library
- **p90 similarity:** top 10% of most-similar pairs (potential duplicates)
- **Histograms:** visualize library diversity and clustering of related designs

> [!NOTE]
> The `jaccard_hist.*` files are only generated if the library includes sequences of differing lengths.

---

## ⚙️ Maintainers
 
- Karl Romanowicz (karl@synplexity.com)
- Calin Plesa (calin@synplexity.com)