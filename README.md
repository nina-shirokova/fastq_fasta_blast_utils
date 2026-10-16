# Bioinformatics Utilities Toolkit

A command-line interface (CLI) tool written in Python for processing bioinformatics file formats (`FASTQ`, `FASTA`, and text-based `BLAST` reports).

## Project Structure

* **`bioinf_utils.py`** — Main execution script containing the CLI router (`argparse`), high-level pipeline orchestration, and directory routing.

* **`helpers/fastq.py`** — Low-level FASTQ parser (streaming via generators) and modular validation functions for GC-content, length, and quality.

* **`helpers/bio_files.py`** — Low-level file I/O validations, FASTA line readers, and text-based BLAST report parsers.

## Features & Capabilities

### 1. FASTQ Filtration (`filter_fastq`)

Filters sequencing reads from a FASTQ file based on customizable biological criteria:

* **GC-content bounds** (`-gc`): Supports single numeric upper limits or custom intervals (e.g., `(20, 80)`). Default value: `(0, 100)`.

* **Length bounds** (`-len`): Filters reads by nucleotide sequence length. Default value: `(0, 2**32)`.

* **Quality threshold** (`-qual`): Discards reads with an average Phred33 quality score below the specified threshold. Default value: `0.0`.

* **Output**: Automatically saves the filtered dataset into the `filtered/` directory.

### 2. Multiline FASTA to One-line Conversion (`convert_fasta`)

* Converts multiline sequence records into a clean single-line format where each sequence occupies exactly one continuous line.

* **Output**: Automatically saves results into the `result/` directory. If the custom output filename (`-o`) is omitted, it automatically generates a name using the prefix `converted_` combined with the input file's base name.

### 3. Text BLAST Output Parsing (`parse_blast`)

* Parses standard text-based BLAST reports.

* Automatically skips multi-line table headers and empty padding to extract the first significant alignment description for each query block.

* Sorts the resulting list of descriptions alphabetically and saves them as a single column.

* **Output**: Automatically saves results into the `result/` directory. If the custom output filename (`-o`) is omitted, it automatically generates a name using the prefix `parsed_` combined with the input file's base name.

## Usage Guide

The toolkit is operated via a unified command-line interface using subcommands.

### General Syntax

```
python bioinf_utils.py <command> [arguments]

```

### Examples

#### 1. Filter FASTQ Reads

```
python bioinf_utils.py filter_fastq -fq input_data.fastq -out filtered_reads.fastq -gc (30, 70) -len (50, 512) -qual 5

#Or using default values:
python bioinf_utils.py filter_fastq -fq input_data.fastq -out filtered_reads.fastq

```

#### 2. Convert Multiline FASTA

```
python bioinf_utils.py convert_fasta -i input_sequences.fasta -o converted_sequences.fasta

# Or using automatic output naming:
python bioinf_utils.py convert_fasta -i input_sequences.fasta

```

#### 3. Parse BLAST Report

```
python bioinf_utils.py parse_blast -i blast_results.txt -o parsed_hits.txt

# Or using automatic output naming:
python bioinf_utils.py parse_blast -i blast_results.txt

```

## Requirements

* Python 3.7 or higher.

* Built entirely using Python standard library modules (`argparse`, `os`) — no external dependencies required.