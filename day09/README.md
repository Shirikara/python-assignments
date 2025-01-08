# Day09 programs
There are two programs that analyze FASTA files and find the longest repeated subsequence that appears twice.

* **analyze_simple.py** requires the path to a file in Fasta or GeneBank format of a file that is physically saved on your PC.
  For example: MIF_FASTA.txt uploaded in this folder.
* **analyze.py** is more sophisticated since it is user-interactive and enables NCBI search, including the new feature called tandem repeats (explained below)

This program allows you to fetch a gene sequence from NCBI, find the longest repeated subsequence, and search for tandem repeats within the gene sequence. 
The program integrates with NCBI's Entrez system to fetch sequences based on a gene name or accession number.

## Features (analyze.py)

- Fetches gene sequences from NCBI using the `gene_name` or `accession number`.
- Identifies the longest repeated subsequence within the sequence.
- Finds and prints tandem repeats of the gene sequence (optional).
- Fully customizable with email address required for NCBI compliance.

## What is Tandem Repeats?
Tandem repeats, also known as Short Tandem Repeats (STRs) or Microsatellites, refer to sequences of DNA where a short motif or pattern of nucleotides is repeated directly adjacent to each other. These repeats can vary in length and number of copies, and they are found at specific loci throughout the genome.

A typical tandem repeat might look like this:

Motif: ATG
Repeated Sequence: ATGATGATG (the motif ATG is repeated three times)
[Tandem Repeats - Wikipedia](https://en.wikipedia.org/wiki/Tandem_repeat)

## Prerequisites

- Python 3.x
- Biopython library (for interacting with NCBI and processing sequences)

You can install the required library using:

```bash
pip install biopython

---
How to Run the Program
Clone the repository or download the analyze.py file.
Open a terminal window and navigate to the directory where the file is saved.
Run the following command:
