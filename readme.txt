# SamReader
Authors : Lucien Maurau, Mathilde Chatain
Contacts :lucien.maurau@etu.umontpellier.fr, mathilde.chatain@etu.umontpellier.fr

# Description

SamReader.py is a Python script designed to analyze SAM (Sequence Alignment/Map) files.
It identifies and categorizes sequencing reads based on their mapping status (unmapped, partially mapped, perfectly mapped), 
analyzes CIGAR strings to compute mutation percentages, and generates multiple summary output files.

# Features
1. SAM file verification
  -check if the file exist
  -ensure the file is not empty
2. Filters reads by mapping quality
3. Identifies:
  -Unmapped reads
  -Partially mapped reads
  -Perfectly mapped reads
4. Handles paired-end reads
5. Parses and summarizes CIGAR strings
6. Computes mutation percentages per read and globally
7. Generates summary and FASTA output files

# Requirements

- Python 3.x
- Standard Python modules : `os`, `sys`, `re`, `itertools`

# Input

A SAM file (`.sam`) produced by `samtools` and a minimum quality value

# Usage

```bash
python3 SamReader.py <input.sam> <min_quality>

# Output Files

All results are written to the `../Results/` directory.

- `only_unmapped.fasta` – unmapped reads
- `only_partially_mapped.fasta` – partially mapped reads
- `summary_unmapped.txt` – count of unmapped reads
- `summary_partially_mapped.txt` – count of partially mapped reads
- Final_Cigar_table.txt - global summary of CIGAR operations and mutation percentages
- QualityDistribution.tsv - distribution of read quality scores