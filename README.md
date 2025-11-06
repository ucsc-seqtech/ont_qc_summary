# Sequencing Quality and Coverage Analysis

## Overview
This project evaluates sequencing data from **Oxford Nanopore (ONT)** focusing on read quality, coverage, and alignment metrics.  
Custom parsing and summary scripts are used to extract coverage statistics, analyze quality bins, and visualize differences between **reported quality** and **alignment-based quality**.

---

## Background
Different sequencing technologies offer tradeoffs between **read length** and **quality**:

| Platform  | Typical Read Length | Quality | Notes |
|------------|--------------------|----------|-------|
| **ONT (Oxford Nanopore)** | 15–20 kb (can reach >100 kb) | Historically lower, improving rapidly | Ideal for long-range structural data |
| **Illumina** | 150–300 bp | Very high | Ideal for short, accurate reads |
| **PacBio** | Up to 100 kb | High | Combines long reads and high accuracy |

---

## Current Setup

### Data Flow
1. **Signal → Basecalling**  
   Raw ONT signal data is converted into basecalled reads (FASTQ/BAM) using [**Dorado**](https://software-docs.nanoporetech.com/dorado/latest/).  

2. **Generate Sequencing Summary**  
   A per-read summary file is created with Dorado using:  `dorado summary <bam> > summary.tsv`

3. **File Parsing**  
   Python scripts (e.g., `calculate_summary_stats.py`) extract:
   - Total read count and total bases  
   - Read length distribution  
   - Quality metrics (mean Q-scores, quality bins)  
   - Total coverage (total Gb / 3.3 Gb ≈ human genome coverage)

4. **Alignment and Quality Evaluation**  
   - Reads are aligned to the reference genome using **minimap2** or **STAR**.  
   - **Alignment quality (mapQ)** is compared to **reported quality (Q-score)** to assess true sequencing performance.

5. **Metadata Tracking**  
   Each dataset includes:
   - Sequencing date  
   - Kit used  
   - Project name  
   - Sample ID  
   - Flowcell ID  

---

## Proposed Work

### 1. Expand Quality Metric Extraction
- Update parsing scripts to include:
  - Coverage bins by mean Q score (e.g., Q20, Q25, Q30)  
  - Average read length and N50 metrics  
  - Mean mapQ for alignment-based quality

### 2. Integrate Alignment-Based Quality
- Compare **reported (per-read)** vs. **true (mapped)** quality metrics.  
- Generate per-sample summaries for both sets of quality measures.

### 3. Generate Summary Tables
- Produce `.txt` summaries with columns for:
  - Read length  
  - Reported quality  
  - Alignment (mapQ) quality  
  - Coverage bins (> Q20, > Q30, etc.)  

### 4. Visualization
- Create plots that display:
  - Coverage vs. read length  
  - Coverage vs. quality bins  
  - Reported vs. alignment quality comparisons  

### 5. Iterative Testing
- Begin with a small subset of BAM files.  
- Scale up once pipeline efficiency is validated.

### 6. Deliverables
- Enhanced `summary_stats_v3.py`  
- Per-sample and aggregate summary tables  
- Visualization scripts (e.g., Python/Matplotlib, R ggplot2)  
- Documentation for sequencing tech center workflows  

## Example Data
- Sequencing summary files generated from aligned bams: [**GoogleDrive**](https://drive.google.com/drive/folders/1vwO2nlkDNYBZ1vg7pxs5SanYoAGv0AdB?usp=sharing).


