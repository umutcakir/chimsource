# ChiMSource Web GUI

A browser-based graphical interface for the ChiMSource pipeline.

## Reference

Çakır, U., Gabed, N., Yurtseven, A., & Kryvoruchko, I. (2025)
*ChiMSource improves the accuracy of studies on novel amino acid sequences
by predicting alternative sources of mass spectrometry-derived peptides.*
Computational and Structural Biotechnology Journal, 27, 3704–3709.
https://doi.org/10.1016/j.csbj.2025.08.023

## Quick Start

```bash
# 1. Install dependencies (one time)
pip install -r requirements.txt

# 2. Launch the server
python app.py

# 3. Open your browser
#    http://localhost:5000
```

The GUI will be accessible from any browser on the same machine.
To change the port, set the environment variable `CHIMSOURCE_PORT`:

```bash
# Windows
set CHIMSOURCE_PORT=8080
python app.py

# macOS / Linux
CHIMSOURCE_PORT=8080 python app.py
```

## Features

- **Drag-and-drop file upload** for nucleotide and peptide FASTA files
- **Real-time progress bars** for file parsing, forward analysis, and reverse complement analysis
- **Live log feed** showing each stage of the pipeline
- **Interactive results table** with pagination and preview of results
- **Download buttons** for the full TSV results and the parameters log
- **All parameters configurable** via the sidebar panel (max gap, codon table, threads, BLAST settings, flanking sequence length, etc.)

## Requirements

- Python 3.9+
- Flask
- pandas
- biotite

## Files

| File | Description |
|------|-------------|
| `app.py` | Flask web server with SSE progress streaming |
| `engine.py` | ChiMSource analysis engine (decoupled from CLI) |
| `templates/index.html` | Web GUI (HTML + CSS + JS, single file) |
| `requirements.txt` | Python dependencies |
