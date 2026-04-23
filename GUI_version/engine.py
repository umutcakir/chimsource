"""
ChiMSource Engine — core analysis logic decoupled from CLI/GUI.

Based on: Çakır, U., Gabed, N., Yurtseven, A., & Kryvoruchko, I. (2025a)
ChiMSource improves the accuracy of studies on novel amino acid sequences by
predicting alternative sources of mass spectrometry-derived peptides.
Computational and Structural Biotechnology Journal, 27, 3704–3709.
"""

import time
import os
import pandas as pd
import re
import multiprocessing
import pickle
import gzip
from functools import partial

from codon_tables import load_codon_table, NCBI_CODON_TABLE_NAMES

# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def find_all_occurrences(translation, amino_acid_sequence):
    indices = []
    index = translation.find(amino_acid_sequence)
    while index != -1:
        indices.append(index)
        index = translation.find(amino_acid_sequence, index + 1)
    return indices


def find_translation_position(nucleotide_sequence, amino_acid_sequence, codon_table):
    positions = []
    for frame in range(3):
        translation = ""
        for i in range(frame, len(nucleotide_sequence) - 2, 3):
            codon = nucleotide_sequence[i:i + 3]
            amino_acid = codon_table.get(codon, "-")
            translation += amino_acid
        if amino_acid_sequence in translation:
            position_list = find_all_occurrences(translation, amino_acid_sequence)
            position = [(index * 3) + frame for index in position_list]
            positions.append(position)
    return [item for sublist in positions for item in sublist]


def translate_sequence(dna_sequence, codon_table):
    protein_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if codon in codon_table:
            protein_sequence += codon_table[codon]
        else:
            protein_sequence += "X"
    return protein_sequence


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))


def clean_header(header):
    return re.sub(r'[\t\n\r]', ' ', header)

# ---------------------------------------------------------------------------
# FASTA parsing  (streaming, with progress callback)
# ---------------------------------------------------------------------------

def _open_fasta(file_path, mode="rb"):
    """Open a FASTA file, transparently handling gzip compression."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def parse_fasta(file_path, progress_callback=None):
    """Parse FASTA (plain or gzipped) with optional progress_callback(bytes_read, total_bytes, n_seqs)."""
    headers, sequences = [], []
    file_size = os.path.getsize(file_path)
    is_gz = file_path.endswith('.gz')
    bytes_read = 0
    current_header = ''
    current_parts = []
    report_interval = max(1, file_size // 200)   # report ~200 times
    last_report = 0

    with _open_fasta(file_path, "rb") as fh:
        for raw_line in fh:
            # Track compressed bytes for gz, raw bytes otherwise
            if is_gz:
                bytes_read = fh.fileobj.tell()
            else:
                bytes_read += len(raw_line)
            line = raw_line.decode('utf-8', errors='replace').strip()
            if line.startswith('>'):
                if current_parts:
                    headers.append(current_header[1:])
                    sequences.append(''.join(current_parts))
                    current_parts = []
                current_header = line
            elif line:
                current_parts.append(line)
            if progress_callback and (bytes_read - last_report >= report_interval):
                progress_callback(bytes_read, file_size, len(headers))
                last_report = bytes_read
        if current_parts:
            sequences.append(''.join(current_parts))
            headers.append(current_header[1:])

    if progress_callback:
        progress_callback(file_size, file_size, len(headers))
    return headers, sequences


def is_valid_fasta(file_path):
    try:
        with _open_fasta(file_path, "rb") as f:
            first_line = f.readline().decode('utf-8', errors='replace').strip()
            return first_line.startswith(">")
    except Exception:
        return False

# ---------------------------------------------------------------------------
# Trimming helpers
# ---------------------------------------------------------------------------

def get_trimmed_sequence(sequence, segment1, segment2, frameshift_pos, seg2_start,
                         max_transcript_length, max_flanking_seq):
    original_length = len(sequence)
    truncated = original_length > max_transcript_length
    if not truncated:
        return sequence, truncated
    segment1_length_nt = len(segment1) * 3
    segment2_length_nt = len(segment2) * 3
    match_start = min(frameshift_pos - segment1_length_nt, seg2_start)
    match_end = max(frameshift_pos, seg2_start + segment2_length_nt)
    start_index = max(0, match_start - max_flanking_seq)
    end_index = min(original_length, match_end + max_flanking_seq)
    return sequence[start_index:end_index], truncated


def find_frameshift_and_segment_starts(sequence, segment1, segment2, offset, codon_table):
    offset1, offset2 = int(offset[0] - 1), int(offset[1] - 1)
    frameshift_happen_list = [
        x + int(len(segment1) * 3) + offset1
        for x in find_translation_position(sequence[offset1:], segment1, codon_table)
    ]
    frameshift_happen_list = [x for x in frameshift_happen_list if (x - offset1) % 3 == 0]
    segment2_start_list = [
        x + offset2
        for x in find_translation_position(sequence[offset2:], segment2, codon_table)
    ]
    segment2_start_list = [x for x in segment2_start_list if (x - offset2) % 3 == 0]
    return frameshift_happen_list, segment2_start_list


def find_frameshift_and_segment_starts_fast(segment1, segment2, offset, frames):
    """Optimized: reuses pre-computed frame translations instead of re-translating."""
    offset1, offset2 = int(offset[0] - 1), int(offset[1] - 1)
    frame1_trans = frames[offset[0]]
    frame2_trans = frames[offset[1]]
    seg1_aa_positions = find_all_occurrences(frame1_trans, segment1)
    frameshift_happen_list = [aa_idx * 3 + len(segment1) * 3 + offset1
                              for aa_idx in seg1_aa_positions]
    seg2_aa_positions = find_all_occurrences(frame2_trans, segment2)
    segment2_start_list = [aa_idx * 3 + offset2
                           for aa_idx in seg2_aa_positions]
    return frameshift_happen_list, segment2_start_list

# ---------------------------------------------------------------------------
# Core analysis  (single nucleotide sequence vs all peptides)
# ---------------------------------------------------------------------------

def process_sequences(fasta1_headers_all, fasta1_sequences_all, fasta2_headers,
                      fasta2_sequences, max_gap,
                      codon_table, max_transcript_length, max_flanking_seq,
                      progress_queue=None, verbose=False, log_queue=None):
    results_all = pd.DataFrame()

    def vlog(msg):
        if verbose and log_queue is not None:
            try:
                log_queue.put(msg)
            except:
                pass

    # Outer loop: iterate over nucleotide sequences
    for b in range(len(fasta1_headers_all)):
        sequence1 = fasta1_sequences_all[b]
        nuc_header = fasta1_headers_all[b]

        # Short name for log messages — extract first word of header
        nuc_short = nuc_header.split()[0][:30] if nuc_header else f"seq_{b+1}"

        # Translate this nucleotide sequence ONCE — reuse for all peptides
        if verbose:
            t_trans = time.time()
            vlog(f"🧬 {nuc_short} ({len(sequence1):,} nt) — translating...")

        frame1 = translate_sequence(sequence1, codon_table)
        frame2 = translate_sequence(sequence1[1:], codon_table)
        frame3 = translate_sequence(sequence1[2:], codon_table)
        frames_dict = {1: frame1, 2: frame2, 3: frame3}

        if verbose:
            vlog(f"🧬 {nuc_short} — translated in {time.time() - t_trans:.1f}s "
                 f"({len(frame1):,} + {len(frame2):,} + {len(frame3):,} aa)")

        # Inner loop: iterate over peptides, reusing the pre-computed frames
        for a in range(len(fasta2_headers)):
            sequence2 = fasta2_sequences[a]

            if verbose:
                pep_t0 = time.time()
                pep_header_short = fasta2_headers[a][:40]
                vlog(f"  {nuc_short} × Peptide {a+1}/{len(fasta2_headers)} ({pep_header_short}, {len(sequence2)} aa)")

            results = {k: [] for k in [
                'Type', 'Frameshift Position', 'Segment 1', 'Segment 2', 'Gap',
                'Frameshift Direction', 'Nucleotide Title', 'Nucleotide Sequence',
                'Protein Title', 'Protein Sequence', 'Truncation for Nucleotide Sequence'
            ]}

            # --- Frameshift search ---
            if verbose:
                t_fs = time.time()
                n_fs_found = 0
            for i in range(1, len(sequence2)):
                seg1 = sequence2[:i]
                seg2 = sequence2[i:]
                if not seg1 or not seg2:
                    continue
                frame_combinations = [
                    (frame1, frame2, [1, 2], 'Frame 1 -> Frame 2'),
                    (frame1, frame3, [1, 3], 'Frame 1 -> Frame 3'),
                    (frame2, frame3, [2, 3], 'Frame 2 -> Frame 3'),
                    (frame2, frame1, [2, 1], 'Frame 2 -> Frame 1'),
                    (frame3, frame1, [3, 1], 'Frame 3 -> Frame 1'),
                    (frame3, frame2, [3, 2], 'Frame 3 -> Frame 2'),
                ]
                for frameA, frameB, offset_val, direction in frame_combinations:
                    if seg1 in frameA and seg2 in frameB:
                        fs_list, s2_list = find_frameshift_and_segment_starts_fast(
                            seg1, seg2, offset_val, frames_dict)
                        for frameshift_pos in fs_list:
                            for seg2_start in s2_list:
                                gap = seg2_start - frameshift_pos
                                if abs(gap) > max_gap:
                                    continue
                                trimmed_seq, truncated = get_trimmed_sequence(
                                    sequence1, seg1, seg2,
                                    frameshift_pos, seg2_start,
                                    max_transcript_length, max_flanking_seq)
                                results['Type'].append('Frameshift')
                                results['Frameshift Position'].append(i)
                                results['Segment 1'].append(seg1)
                                results['Segment 2'].append(seg2)
                                results['Gap'].append(f"{gap:+d}")
                                results['Frameshift Direction'].append(direction)
                                results['Nucleotide Title'].append(nuc_header)
                                results['Nucleotide Sequence'].append(trimmed_seq)
                                results['Protein Title'].append(fasta2_headers[a])
                                results['Protein Sequence'].append(sequence2)
                                results['Truncation for Nucleotide Sequence'].append(str(truncated))
                                if verbose:
                                    n_fs_found += 1

            if verbose:
                vlog(f"    {nuc_short} × P{a+1}: Frameshift {time.time() - t_fs:.1f}s  "
                     f"({len(sequence2)-1} splits × 6 combos, {n_fs_found} hits)")

            # --- Direct (no frameshift) search  [FIXED: match-centred trimming] ---
            if verbose:
                t_dm = time.time()
                n_dm_found = 0
            frames = {
                'frame1': (frame1, 0),
                'frame2': (frame2, 1),
                'frame3': (frame3, 2),
            }
            for frame_name, (frame_content, frame_offset) in frames.items():
                if sequence2 in frame_content:
                    aa_positions = find_all_occurrences(frame_content, sequence2)
                    for aa_index in aa_positions:
                        truncated = len(sequence1) > max_transcript_length
                        if truncated:
                            match_start_nt = aa_index * 3 + frame_offset
                            match_end_nt = match_start_nt + len(sequence2) * 3
                            start_idx = max(0, match_start_nt - max_flanking_seq)
                            end_idx = min(len(sequence1), match_end_nt + max_flanking_seq)
                            trimmed_seq = sequence1[start_idx:end_idx]
                        else:
                            trimmed_seq = sequence1
                        results['Type'].append('Without frameshift')
                        results['Frameshift Position'].append("No Frameshift")
                        results['Segment 1'].append(sequence2)
                        results['Segment 2'].append(" ")
                        results['Gap'].append("+0")
                        results['Frameshift Direction'].append(
                            f'The sequence is present in {frame_name} without any frameshift')
                        results['Nucleotide Title'].append(nuc_header)
                        results['Nucleotide Sequence'].append(trimmed_seq)
                        results['Protein Title'].append(fasta2_headers[a])
                        results['Protein Sequence'].append(sequence2)
                        results['Truncation for Nucleotide Sequence'].append(str(truncated))
                        if verbose:
                            n_dm_found += 1

            if verbose:
                vlog(f"    {nuc_short} × P{a+1}: Direct {time.time() - t_dm:.1f}s ({n_dm_found} hits)")
                vlog(f"    {nuc_short} × P{a+1}: Done in {time.time() - pep_t0:.1f}s")

            results_df = pd.DataFrame(results)
            results_all = pd.concat([results_all, results_df], ignore_index=True)

            # Signal per-peptide progress
            if progress_queue is not None:
                try:
                    progress_queue.put('peptide_done')
                except:
                    pass

    return results_all

# ---------------------------------------------------------------------------
# Parallel wrapper — chunk-based, no Manager
# ---------------------------------------------------------------------------

_shared_gui = {}

def _gui_worker_init(pep_headers, pep_sequences, codon_table, max_gap,
                     max_transcript_length, max_flanking_seq,
                     verbose, progress_counter, progress_lock,
                     log_queue, input_profile):
    """Called once per worker — store shared read-only data."""
    _shared_gui['pep_headers'] = pep_headers
    _shared_gui['pep_sequences'] = pep_sequences
    _shared_gui['codon_table'] = codon_table
    _shared_gui['max_gap'] = max_gap
    _shared_gui['max_transcript_length'] = max_transcript_length
    _shared_gui['max_flanking_seq'] = max_flanking_seq
    _shared_gui['verbose'] = verbose
    _shared_gui['progress_counter'] = progress_counter
    _shared_gui['progress_lock'] = progress_lock
    _shared_gui['log_queue'] = log_queue
    _shared_gui['input_profile'] = input_profile


def _gui_process_chunk(chunk):
    """Process a list of (header, sequence) pairs — one chunk per worker."""
    pep_headers = _shared_gui['pep_headers']
    pep_sequences = _shared_gui['pep_sequences']
    codon_table = _shared_gui['codon_table']
    max_gap = _shared_gui['max_gap']
    max_transcript_length = _shared_gui['max_transcript_length']
    max_flanking_seq = _shared_gui['max_flanking_seq']
    verbose = _shared_gui['verbose']
    counter = _shared_gui['progress_counter']
    lock = _shared_gui['progress_lock']
    log_q = _shared_gui['log_queue']
    profile = _shared_gui['input_profile']

    genomic_verbose = verbose and profile == 'genomic'
    trans_verbose = verbose and profile == 'transcriptomic'

    def vlog(msg):
        if log_q is not None:
            try:
                log_q.put(msg)
            except:
                pass

    all_results = []
    local_progress = 0
    FLUSH_EVERY = 50
    trans_seq_count = 0
    trans_hit_count = 0
    trans_chunk_start = time.time()

    for nuc_header, sequence1 in chunk:
        nuc_short = nuc_header.split()[0][:30] if nuc_header else "seq"

        if genomic_verbose:
            t_trans = time.time()
            vlog(f"🧬 {nuc_short} ({len(sequence1):,} nt) — translating...")

        frame1 = translate_sequence(sequence1, codon_table)
        frame2 = translate_sequence(sequence1[1:], codon_table)
        frame3 = translate_sequence(sequence1[2:], codon_table)
        frames_dict = {1: frame1, 2: frame2, 3: frame3}

        if genomic_verbose:
            vlog(f"🧬 {nuc_short} — translated in {time.time() - t_trans:.1f}s "
                 f"({len(frame1):,} + {len(frame2):,} + {len(frame3):,} aa)")

        for a in range(len(pep_headers)):
            sequence2 = pep_sequences[a]

            if genomic_verbose:
                pep_t0 = time.time()
                pep_short = pep_headers[a][:40]
                vlog(f"  {nuc_short} × Peptide {a+1}/{len(pep_headers)} ({pep_short}, {len(sequence2)} aa)")

            results = {k: [] for k in [
                'Type', 'Frameshift Position', 'Segment 1', 'Segment 2', 'Gap',
                'Frameshift Direction', 'Nucleotide Title', 'Nucleotide Sequence',
                'Protein Title', 'Protein Sequence', 'Truncation for Nucleotide Sequence'
            ]}

            # --- Frameshift search ---
            if genomic_verbose:
                t_fs = time.time()
                n_fs_found = 0
            for i in range(1, len(sequence2)):
                seg1 = sequence2[:i]
                seg2 = sequence2[i:]
                if not seg1 or not seg2:
                    continue
                frame_combinations = [
                    (frame1, frame2, [1, 2], 'Frame 1 -> Frame 2'),
                    (frame1, frame3, [1, 3], 'Frame 1 -> Frame 3'),
                    (frame2, frame3, [2, 3], 'Frame 2 -> Frame 3'),
                    (frame2, frame1, [2, 1], 'Frame 2 -> Frame 1'),
                    (frame3, frame1, [3, 1], 'Frame 3 -> Frame 1'),
                    (frame3, frame2, [3, 2], 'Frame 3 -> Frame 2'),
                ]
                for frameA, frameB, offset_val, direction in frame_combinations:
                    if seg1 in frameA and seg2 in frameB:
                        fs_list, s2_list = find_frameshift_and_segment_starts_fast(
                            seg1, seg2, offset_val, frames_dict)
                        for frameshift_pos in fs_list:
                            for seg2_start in s2_list:
                                gap = seg2_start - frameshift_pos
                                if abs(gap) > max_gap:
                                    continue
                                trimmed_seq, truncated = get_trimmed_sequence(
                                    sequence1, seg1, seg2,
                                    frameshift_pos, seg2_start,
                                    max_transcript_length, max_flanking_seq)
                                results['Type'].append('Frameshift')
                                results['Frameshift Position'].append(i)
                                results['Segment 1'].append(seg1)
                                results['Segment 2'].append(seg2)
                                results['Gap'].append(f"{gap:+d}")
                                results['Frameshift Direction'].append(direction)
                                results['Nucleotide Title'].append(nuc_header)
                                results['Nucleotide Sequence'].append(trimmed_seq)
                                results['Protein Title'].append(pep_headers[a])
                                results['Protein Sequence'].append(sequence2)
                                results['Truncation for Nucleotide Sequence'].append(str(truncated))
                                if genomic_verbose:
                                    n_fs_found += 1

            if genomic_verbose:
                vlog(f"    {nuc_short} × P{a+1}: Frameshift {time.time() - t_fs:.1f}s  "
                     f"({len(sequence2)-1} splits × 6 combos, {n_fs_found} hits)")

            # --- Direct match search ---
            if genomic_verbose:
                t_dm = time.time()
                n_dm_found = 0
            frames = {
                'frame1': (frame1, 0),
                'frame2': (frame2, 1),
                'frame3': (frame3, 2),
            }
            for frame_name, (frame_content, frame_offset) in frames.items():
                if sequence2 in frame_content:
                    aa_positions = find_all_occurrences(frame_content, sequence2)
                    for aa_index in aa_positions:
                        truncated = len(sequence1) > max_transcript_length
                        if truncated:
                            match_start_nt = aa_index * 3 + frame_offset
                            match_end_nt = match_start_nt + len(sequence2) * 3
                            start_idx = max(0, match_start_nt - max_flanking_seq)
                            end_idx = min(len(sequence1), match_end_nt + max_flanking_seq)
                            trimmed_seq = sequence1[start_idx:end_idx]
                        else:
                            trimmed_seq = sequence1
                        results['Type'].append('Without frameshift')
                        results['Frameshift Position'].append("No Frameshift")
                        results['Segment 1'].append(sequence2)
                        results['Segment 2'].append(" ")
                        results['Gap'].append("+0")
                        results['Frameshift Direction'].append(
                            f'The sequence is present in {frame_name} without any frameshift')
                        results['Nucleotide Title'].append(nuc_header)
                        results['Nucleotide Sequence'].append(trimmed_seq)
                        results['Protein Title'].append(pep_headers[a])
                        results['Protein Sequence'].append(sequence2)
                        results['Truncation for Nucleotide Sequence'].append(str(truncated))
                        if genomic_verbose:
                            n_dm_found += 1

            if genomic_verbose:
                vlog(f"    {nuc_short} × P{a+1}: Direct {time.time() - t_dm:.1f}s ({n_dm_found} hits)")
                vlog(f"    {nuc_short} × P{a+1}: Done in {time.time() - pep_t0:.1f}s")

            n_hits_this_pair = len(results['Type'])
            if n_hits_this_pair:
                all_results.append(pd.DataFrame(results))
                if trans_verbose:
                    trans_hit_count += n_hits_this_pair

            local_progress += 1

            if local_progress >= FLUSH_EVERY:
                with lock:
                    counter.value += local_progress
                local_progress = 0

        # End-of-transcript flush
        if local_progress > 0:
            with lock:
                counter.value += local_progress
            local_progress = 0

        # Transcriptomic verbose: periodic summary every 500 sequences
        if trans_verbose:
            trans_seq_count += 1
            if trans_seq_count % 500 == 0:
                elapsed = time.time() - trans_chunk_start
                rate = trans_seq_count / max(elapsed, 0.001)
                remaining = len(chunk) - trans_seq_count
                eta = remaining / max(rate, 0.001)
                vlog(f"📊 Worker: {trans_seq_count}/{len(chunk)} transcripts "
                     f"({rate:.1f} seq/s, ~{eta/60:.1f} min remaining, "
                     f"{trans_hit_count} hits)")

    # Transcriptomic verbose: final summary
    if trans_verbose:
        elapsed = time.time() - trans_chunk_start
        vlog(f"✅ Worker done: {trans_seq_count} transcripts in {elapsed/60:.1f} min "
             f"({trans_seq_count/max(elapsed,0.001):.1f} seq/s, {trans_hit_count} hits)")

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


def process_sequences_parallel(fasta1_headers, fasta1_sequences,
                               fasta2_headers, fasta2_sequences,
                               max_gap,
                               n_jobs, codon_table,
                               max_transcript_length, max_flanking_seq,
                               progress_callback=None, verbose=False,
                               log_callback=None):
    """Chunk-based parallelization with adaptive presentation."""
    import pickle as _pickle

    n_nuc = len(fasta1_headers)
    n_pep = len(fasta2_headers)
    total_pairs = n_nuc * n_pep

    # Detect input profile
    avg_seq_len = sum(len(s) for s in fasta1_sequences) / max(n_nuc, 1)
    if n_nuc <= 200 and avg_seq_len > 500_000:
        input_profile = 'genomic'
    else:
        input_profile = 'transcriptomic'

    profile_label = "Genomic" if input_profile == 'genomic' else "Transcriptomic"
    if log_callback:
        log_callback(f"📊 Input profile: {profile_label} "
                     f"({n_nuc:,} sequences, avg {avg_seq_len/1000:.1f} KB × {n_pep} peptides = "
                     f"{total_pairs:,} pairs)")

    progress_counter = multiprocessing.Value('i', 0)
    progress_lock = multiprocessing.Lock()

    # Only create a Manager for verbose log queue (lightweight — few messages)
    mgr = None
    log_queue = None
    if verbose and log_callback:
        mgr = multiprocessing.Manager()
        log_queue = mgr.Queue()

    nuc_pairs = list(zip(fasta1_headers, fasta1_sequences))
    actual_workers = min(n_jobs, n_nuc)
    chunk_size = max(1, (n_nuc + actual_workers - 1) // actual_workers)
    chunks = [nuc_pairs[i:i + chunk_size] for i in range(0, n_nuc, chunk_size)]

    pairs_done = 0

    with multiprocessing.Pool(
        processes=actual_workers,
        initializer=_gui_worker_init,
        initargs=(fasta2_headers, fasta2_sequences, codon_table, max_gap,
                  max_transcript_length, max_flanking_seq,
                  verbose, progress_counter, progress_lock,
                  log_queue, input_profile)
    ) as pool:
        async_result = pool.map_async(_gui_process_chunk, chunks)

        while not async_result.ready():
            # Drain log queue
            if log_queue is not None:
                while not log_queue.empty():
                    try:
                        msg = log_queue.get_nowait()
                        if log_callback:
                            log_callback(msg)
                    except:
                        break
            # Poll progress counter
            time.sleep(0.3)
            current = progress_counter.value
            if current > pairs_done:
                if progress_callback:
                    progress_callback(current, total_pairs)
                pairs_done = current

        # Final drain
        if log_queue is not None:
            while not log_queue.empty():
                try:
                    msg = log_queue.get_nowait()
                    if log_callback:
                        log_callback(msg)
                except:
                    break
        current = progress_counter.value
        if current > pairs_done and progress_callback:
            progress_callback(current, total_pairs)

        chunk_results = async_result.get()

    non_empty = [r for r in chunk_results if r is not None and not r.empty]
    if non_empty:
        return pd.concat(non_empty, ignore_index=True)
    return pd.DataFrame()

# ---------------------------------------------------------------------------
# High-level run function (used by both CLI and web GUI)
# ---------------------------------------------------------------------------

def run_analysis(nucleotide_file, peptide_file, *,
                 codon_table_id=1, max_gap=2, num_threads=1,
                 output_base='results',
                 max_transcript_length=30000, max_flanking_seq=500,
                 no_reverse_complement_check=False,
                 verbose=False,
                 log_callback=None, progress_callback=None):
    """
    Run the full ChiMSource pipeline.

    log_callback(message: str)         — status messages
    progress_callback(stage, done, total) — numeric progress
        stage is one of: 'parse_nuc', 'parse_pep', 'forward', 'reverse'

    Returns (combined_df, params_text, output_tsv_path)
    """

    def log(msg):
        if log_callback:
            log_callback(msg)

    def prog(stage):
        def _inner(done, total, *extra):
            if progress_callback:
                progress_callback(stage, done, total)
        return _inner

    from datetime import datetime
    analysis_started_at = datetime.now()
    start_time = time.time()

    # --- Codon table ---
    codon_table, ct_name = load_codon_table(codon_table_id)
    log(f"Using codon table: {ct_name}")

    # --- Validate ---
    for path in [nucleotide_file, peptide_file]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")
        if not is_valid_fasta(path):
            raise ValueError(f"Invalid FASTA format: {path}")

    # --- Parse nucleotides ---
    log("Parsing nucleotide file...")
    raw_nuc_h, raw_nuc_s = parse_fasta(nucleotide_file, progress_callback=prog('parse_nuc'))
    nuc_headers = [clean_header(h) for h in raw_nuc_h]
    nuc_seqs = list(raw_nuc_s)
    log(f"Loaded {len(nuc_headers)} nucleotide sequences")

    # --- Parse peptides ---
    log("Parsing peptide file...")
    raw_pep_h, raw_pep_s = parse_fasta(peptide_file, progress_callback=prog('parse_pep'))
    pep_headers = [clean_header(h) for h in raw_pep_h]
    pep_seqs = list(raw_pep_s)
    log(f"Loaded {len(pep_headers)} peptide sequences")

    # --- Forward ---
    log("▶ STAGE 1/2: Forward strand analysis")
    fwd_start = time.time()
    forward_results = process_sequences_parallel(
        nuc_headers, nuc_seqs, pep_headers, pep_seqs,
        max_gap, num_threads, codon_table,
        max_transcript_length, max_flanking_seq,
        progress_callback=prog('forward'), verbose=verbose,
        log_callback=log)
    fwd_elapsed = time.time() - fwd_start
    if not forward_results.empty:
        forward_results['Frame Direction'] = 'Forward'
    log(f"✅ Forward analysis completed in {fwd_elapsed/60:.1f} minutes")

    # --- Reverse complement ---
    reverse_results = pd.DataFrame()
    rev_elapsed = 0.0
    if not no_reverse_complement_check:
        log("◀ STAGE 2/2: Reverse complement analysis")
        rev_start = time.time()
        rev_seqs = [reverse_complement(s) for s in nuc_seqs]
        reverse_results = process_sequences_parallel(
            nuc_headers, rev_seqs, pep_headers, pep_seqs,
            max_gap, num_threads, codon_table,
            max_transcript_length, max_flanking_seq,
            progress_callback=prog('reverse'), verbose=verbose,
            log_callback=log)
        rev_elapsed = time.time() - rev_start
        if not reverse_results.empty:
            reverse_results['Frame Direction'] = 'Reverse'
        log(f"✅ Reverse complement analysis completed in {rev_elapsed/60:.1f} minutes")

    # --- Combine ---
    combined = pd.concat([forward_results, reverse_results], ignore_index=True)
    if not combined.empty:
        # Sort results by peptide input order
        peptide_order = {title: idx for idx, title in enumerate(pep_headers)}
        combined['_pep_order'] = combined['Protein Title'].map(peptide_order)
        combined = combined.sort_values('_pep_order', kind='stable').drop(columns=['_pep_order']).reset_index(drop=True)

        column_order = [
            'Protein Title', 'Protein Sequence', 'Type',
            'Frameshift Position', 'Segment 1', 'Segment 2', 'Gap',
            'Frameshift Direction', 'Nucleotide Title', 'Nucleotide Sequence',
            'Frame Direction', 'Truncation for Nucleotide Sequence'
        ]
        combined = combined[column_order]

    analysis_finished_at = datetime.now()
    elapsed = time.time() - start_time
    elapsed_dt = (analysis_finished_at - analysis_started_at).total_seconds()

    def format_time(secs):
        h = int(secs // 3600)
        m = int((secs % 3600) // 60)
        s = int(secs % 60)
        if h > 0:
            return f"{h}h {m}m {s}s ({secs:.0f} seconds)"
        elif m > 0:
            return f"{m}m {s}s ({secs:.0f} seconds)"
        else:
            return f"{s}s ({secs:.2f} seconds)"

    elapsed_str = format_time(elapsed_dt)

    # --- Save ---
    output_tsv = f"{output_base}.tsv"
    if not combined.empty:
        combined.to_csv(output_tsv, sep='\t', index=False)
    log(f"Analysis complete in {elapsed_str} — {len(combined)} results")

    fwd_str = f"{fwd_elapsed/60:.1f} minutes ({fwd_elapsed:.0f} seconds)"
    rev_str = f"{rev_elapsed/60:.1f} minutes ({rev_elapsed:.0f} seconds)" if rev_elapsed > 0 else "Skipped"

    params_text = f"""Analysis Parameters
{'=' * 18}
Nucleotide file: {nucleotide_file}
Peptide file: {peptide_file}
Codon table: {ct_name}
Max gap size: {max_gap}
Max transcript length: {max_transcript_length}
Max flanking sequence: {max_flanking_seq}
Reverse complement check: {'Disabled' if no_reverse_complement_check else 'Enabled'}
Verbose logging: {'Enabled' if verbose else 'Disabled'}
Number of threads: {num_threads}
Total sequences processed: {len(nuc_headers)} nucleotide sequences, {len(pep_headers)} peptides
Analysis started at: {analysis_started_at.strftime('%Y-%m-%d %H:%M:%S')}
Analysis finished at: {analysis_finished_at.strftime('%Y-%m-%d %H:%M:%S')}
Forward strand analysis time: {fwd_str}
Reverse complement analysis time: {rev_str}
Total analysis time: {elapsed_str}
"""
    params_path = f"{output_base}_parameters.txt"
    with open(params_path, 'w') as f:
        f.write(params_text)

    return combined, params_text, output_tsv
