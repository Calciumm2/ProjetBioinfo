SamReader

Description

SamReader is a Python script designed to analyze SAM (Sequence Alignment/Map) files.
It identifies and categorizes sequencing reads based on their mapping status (unmapped, partially mapped, perfectly mapped), analyzes CIGAR strings to compute mutation percentages, and generates multiple summary output files.

Authors

Lucien Maurau – lucien.maurau@etu.umontpellier.fr
Mathilde Chatain – mathilde.chatain@etu.umontpellier.fr


Requirements

Python 3.x
Standard Python libraries only:

  * `os`
  * `sys`
  * `re`
  * `itertools`

No external dependencies are required.

Input

A SAM file (`.sam`) produced by `samtools` and a minimum of quality

Usage

bash
python3 SamReader.py <input.sam>
```

Output Files

All results are written to the `../Results/` directory.

- `only_unmapped.fasta` – unmapped reads
- `only_partially_mapped.fasta` – partially mapped reads
- `summary_unmapped.txt` – count of unmapped reads
- `summary_partially_mapped.txt` – count of partially mapped reads
- Final_Cigar_table.txt
- QualityDistribution.tsv

Main Features

- Checks input file existence and validity
- Filters reads by mapping quality
- Identifies:
  * Unmapped reads
  * Partially mapped reads
  * Perfectly mapped reads
- Handles paired-end reads
- Parses and summarizes CIGAR strings
- Computes mutation percentages per read and globally
- Generates summary and FASTA output files