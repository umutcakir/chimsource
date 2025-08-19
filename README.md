# ChiMSource: Chimeric MS Peptide Sources

A computational pipeline for detecting alternative chimeric and non-chimeric sources of mass spectrometry-identified peptides.

---

## Overview

ChiMSource enables systematic identification of alternative genomic sources for chimeric peptides. The tool addresses critical challenges in programmed ribosomal frameshifting (PRF) studies by: 

- Identifying potential alternative loci for chimeric peptides
- Differentiating between unique and multi-source chimeric origins
- Detecting reverse complement translation events
- Analyzing repeat element contributions to chimeric translation

Key applications include functional genomics studies of non-canonical translation events and validation of putative PRF sites.

---

## Features


- Multi-frame translation analysis (forward and reverse frames 0, +1, and +2)
- Repeat element integration for comprehensive source detection
- CSV output with frameshift position annotations
- Parallel processing support for large genomes

---

## Available Commands

ChiMSource operates through a single command with multiple configurable parameters:

### chimsource

Identifies chimeric sequences in nucleotide databases through frameshift analysis.

### Usage

> chimsource --nucleotide_file <genome.fasta> --peptide_file <peptides.fasta> [--max_gap 2] [--codon_table 1] [--num_threads 1] [--blast False] [--sensitivity 5.0] [--word_size 2] [--evalue 1000] [--output results] [--max_transcript_length 30000] [--max_flanking_seq 500] [--no_reverse_complement_check False]

### Parameters

- --nucleotide_file (required): Path to nucleotide FASTA file (genome/transcriptome/repeatome).
- --peptide_file: Path to query peptide FASTA file.
- --max_gap: Maximum allowed gap size in nucleotides (default: 2). This parameter corresponds to the frameshift value. Default value of 2 will permit consideration of frameshifts by zero (no frameshift), one, and two nucleotides in either direction. The maximal setting of this parameter is user-defined such as 10, which means unconventionally “long” frameshifting events can be considered as alternative sources of a given query peptide sequence.
- --codon_table: ID of the NCBI codon table (default: The Standard Code). Users can import other codon tables by typing corresponding table IDs (NCBI). The ID of the NCBI codon table can be accessed via the following link: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- --num_threads: Number of parallel processing threads (default: 1). This parameter can be user-defined based on the availability of processors in the system.
- --blast: Enable BLAST pre-filtering for acceleration (not recommended). Omitting this flag will disable BLAST. BLAST pre-filtering can be of use only if the input file is very large and only if losing the information on some alternative sources is acceptable for a specific study. This option is available if BLAST was installed and available in the PATH. If one’s system did not add BLAST to the PATH environment variable, only one folder needs to be added to the system’s PATH, the one containing all the BLAST executables (blastn.exe, tblastn.exe, blastx.exe, etc.).
- --threshold: Word inclusion threshold (higher = slower) for the BLAST acceleration (default: 5.0).
- --word_size: BLAST word size parameter (default: 2).
- --evalue: BLAST e-value threshold (default: 1000).
- --output: Base name for the CSV output file (default: results).
- --max_transcript_length: Max length of a transcript sequence (default: 30,000, in nucleotides). If a subject sequence exceeds this threshold, it will be trimmed according to the setting of [--max_flanking_seq] below.
- --max_flanking_seq: Max length of flanking sequences to include (default: 500, in nucleotides). If a subject sequence is shorter than the value set in [--max_transcript_length] (see above), it will be recorded as a whole sequence in the output file. For subject sequences longer than [--max_transcript_length] (e.g., a chromosome), this setting will record the frameshift region (the complete matching part) plus [--max_flanking_seq] nucleotides upstream and downstream of the matching site.
- --no_reverse_complement_check: Omitting this flag will enable reverse complement analysis. This is recommended because the software does not discriminate between transcript, repeat, and genomic sequences. In genomic and repeat sequences, both strands can encode peptides and proteins, even in the same region (oppositely overlapping genes).

### Output
Contains tabular data on detected chimeric and non-chimeric sequences with these columns:

- Type: Detection type. “Frameshift” stands for a chimeric alternative source. “Without frameshift” stands for a non-chimeric alternative source (direct match).
- Frameshift Position: Position of the last amino acid before the frameshift (position 1 is the first amino acid in a query peptide).
- Segment 1: Contains a portion of a query peptide sequence before the frameshift.
- Segment 2: Contains a portion of a query peptide sequence after the frameshift.
- Gap: Indicates the frameshifts value in the alternative source. If the maximal gap value is set to 2 (default), the following alternative sources will be reported: 0 (no frameshift), +1 and +2 (forward frameshifts), -1 and -2 (backward frameshifts).
- Frameshift Direction: Reading frame transition (Frame n -> Frame m).
- Nucleotide Title: Contains a header of a subject sequence from the input nucleotide FASTA file.
- Nucleotide Sequence: Contains a subject sequence delimited by the setting in [--max_flanking_seq].
- Protein Title: Contains a query header from the input peptide FASTA file.
- Protein Sequence: Contains a query peptide sequence.
- Frame Direction: Strand orientation. Forward frames alone (e.g., transcriptome input) or together with reverse frames (genomic DNA or repeatome input) can be considered.
- Truncation for Nucleotide Sequence (True or False): If the subject nucleotide sequence exceeds [--max_transcript_length], it will be truncated in the given [--max_flanking_seq] threshold level (output True). Else, the output is False.

---

## ChiMSource Installation Guide

### Requirements

- python 3.8+
- biopython >= 1.81
- pandas >= 2.0
- NCBI BLAST+ (tblastn)


### Installation

- Install from PyPI
> pip install chimsource


- For development/editable installation
> git clone https://github.com/umutcakir/chimsource

> cd chimsource

> pip install -e .

ChiMSource was developed to advance research on mosaic translation and programmed ribosomal frameshifting. It has enabled the discrimination between chimeric proteins with unique sources in the genome and proteins that have multiple alternative sources, which is essential for their functional analysis. For biological context and related studies, see our publication: Çakır et al. (2024, preprint).

If you use ChiMSource in your research, please cite the following article:

Umut Çakır, Noujoud Gabed, Ali Yurtseven, Igor Kryvoruchko (2025). ChiMSource, a program that predicts potential alternative sources of mass spectrometry-derived peptides, improves the accuracy of studies on novel amino acid sequences. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2025.05.30.656965

To report bugs, ask questions, or suggest features, feel free to open an issue on GitHub. Your feedback and citations help us improve and sustain this tool.
