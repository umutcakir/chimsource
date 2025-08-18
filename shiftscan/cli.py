import time
import os
import biotite
import biotite.sequence
import pandas as pd
import subprocess
import argparse
from joblib import Parallel, delayed
from tqdm import tqdm
import sys
import re
import pickle
import multiprocessing

# Default codon table
codon_table_default = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


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
            amino_acid = codon_table[codon]
            protein_sequence += amino_acid
        else:
            protein_sequence += "X"  # Unknown codon
    return protein_sequence


def parse_fasta(file_path):
    headers = []
    sequences = []
    with open(file_path, "r") as file:
        lines = file.readlines()
        current_header = ''
        current_sequence = ''
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    headers.append(current_header[1:])  # Remove '>'
                    sequences.append(current_sequence)
                    current_sequence = ''
                current_header = line
            else:
                current_sequence += line
        if current_sequence:  # Capture the last sequence in the file
            sequences.append(current_sequence)
            headers.append(current_header[1:])  # Remove '>'
    return headers, sequences


def get_trimmed_sequence(sequence, segment1, segment2, frameshift_pos, seg2_start,
                         max_transcript_length, max_flanking_seq):
    """Trim sequence to match region with flanking sequences and return truncation status"""
    original_length = len(sequence)
    truncated = original_length > max_transcript_length

    if not truncated:
        return sequence, truncated

    # Calculate nucleotide positions
    segment1_length_nt = len(segment1) * 3
    segment2_length_nt = len(segment2) * 3

    # Calculate start and end of match region
    match_start = min(frameshift_pos - segment1_length_nt, seg2_start)
    match_end = max(frameshift_pos, seg2_start + segment2_length_nt)

    # Expand with flanking sequences
    start_index = max(0, match_start - max_flanking_seq)
    end_index = min(original_length, match_end + max_flanking_seq)

    return sequence[start_index:end_index], truncated


def find_frameshift_and_segment_starts(sequence, segment1, segment2, offset, codon_table):
    offset1, offset2 = int(offset[0] - 1), int(offset[1] - 1)

    frameshift_happen_list = [x + int((len(segment1) * 3)) + offset1
                              for x in find_translation_position(sequence[offset1:], segment1, codon_table)]
    frameshift_happen_list = [x for x in frameshift_happen_list if (x - offset1) % 3 == 0]

    segment2_start_list = [x + offset2
                           for x in find_translation_position(sequence[offset2:], segment2, codon_table)]
    segment2_start_list = [x for x in segment2_start_list if (x - offset2) % 3 == 0]

    return frameshift_happen_list, segment2_start_list


def process_sequences(fasta1_headers_all, fasta1_sequences_all, fasta2_headers,
                      fasta2_sequences, max_gap, speed_up_by_blast, blast_output,
                      codon_table, max_transcript_length, max_flanking_seq):
    results_all = pd.DataFrame()

    for a in range(len(fasta2_headers)):
        if speed_up_by_blast and not blast_output.empty:
            desired_rows = blast_output[blast_output[0] == fasta2_headers[a]]
            unique_nt_entries = desired_rows[1].tolist()
            selected = list(set(unique_nt_entries))

            fasta1_headers = []
            fasta1_sequences = []
            for header, sequence in zip(fasta1_headers_all, fasta1_sequences_all):
                if any(sel in header for sel in selected):
                    fasta1_headers.append(header)
                    fasta1_sequences.append(sequence)
        else:
            fasta1_headers, fasta1_sequences = fasta1_headers_all, fasta1_sequences_all

        for b in range(len(fasta1_headers)):
            sequence2 = fasta2_sequences[a]
            sequence1 = fasta1_sequences[b]

            # Translate sequence1 into all three reading frames
            frame1 = translate_sequence(sequence1, codon_table)
            frame2 = translate_sequence(sequence1[1:], codon_table)
            frame3 = translate_sequence(sequence1[2:], codon_table)

            # Store results
            results = {
                'Type': [],
                'Frameshift Position': [],
                'Segment 1': [],
                'Segment 2': [],
                'Gap': [],
                'Frameshift Direction': [],
                'Nucleotide Title': [],
                'Nucleotide Sequence': [],
                'Protein Title': [],
                'Protein Sequence': [],
                'Truncation for Nucleotide Sequence': []
            }

            # Check for frameshift positions
            for i in range(1, len(sequence2)):  # Start from position 1 (1-based)
                segment1 = sequence2[:i]
                segment2 = sequence2[i:]

                if len(segment1) == 0 or len(segment2) == 0:
                    continue

                # Check all possible frame combinations
                frame_combinations = [
                    (frame1, frame2, [1, 2], 'Frame 1 -> Frame 2'),
                    (frame1, frame3, [1, 3], 'Frame 1 -> Frame 3'),
                    (frame2, frame3, [2, 3], 'Frame 2 -> Frame 3'),
                    (frame2, frame1, [2, 1], 'Frame 2 -> Frame 1'),
                    (frame3, frame1, [3, 1], 'Frame 3 -> Frame 1'),
                    (frame3, frame2, [3, 2], 'Frame 3 -> Frame 2')
                ]

                for frameA, frameB, offset_val, direction in frame_combinations:
                    if segment1 in frameA and segment2 in frameB:
                        frameshift_list, segment2_list = find_frameshift_and_segment_starts(
                            sequence1, segment1, segment2, offset_val, codon_table
                        )
                        for frameshift_pos in frameshift_list:
                            for seg2_start in segment2_list:
                                gap = seg2_start - frameshift_pos
                                if abs(gap) > max_gap:
                                    continue

                                # Get trimmed sequence and truncation status
                                trimmed_seq, truncated = get_trimmed_sequence(
                                    sequence1, segment1, segment2,
                                    frameshift_pos, seg2_start,
                                    max_transcript_length, max_flanking_seq
                                )

                                # Format gap with sign
                                gap_str = f"{gap:+d}"

                                # Add results
                                results['Type'].append('Frameshift')
                                results['Frameshift Position'].append(i)  # 1-based position
                                results['Segment 1'].append(segment1)
                                results['Segment 2'].append(segment2)
                                results['Gap'].append(gap_str)
                                results['Frameshift Direction'].append(direction)
                                results['Nucleotide Title'].append(fasta1_headers[b])
                                results['Nucleotide Sequence'].append(trimmed_seq)
                                results['Protein Title'].append(fasta2_headers[a])
                                results['Protein Sequence'].append(sequence2)
                                results['Truncation for Nucleotide Sequence'].append(str(truncated))

            # Check for direct presence without frameshift
            frames = {
                'frame1': frame1,
                'frame2': frame2,
                'frame3': frame3
            }
            for frame_name, frame_content in frames.items():
                if sequence2 in frame_content:
                    # For direct matches, determine truncation status
                    truncated = len(sequence1) > max_transcript_length

                    # Format gap as 0 with sign
                    gap_str = "+0"

                    # Get sequence to display
                    if truncated:
                        # For long sequences, take a central portion
                        start_idx = max(0, len(sequence1) // 2 - max_flanking_seq)
                        end_idx = min(len(sequence1), len(sequence1) // 2 + max_flanking_seq)
                        trimmed_seq = sequence1[start_idx:end_idx]
                    else:
                        trimmed_seq = sequence1

                    results['Type'].append('Without frameshift')
                    results['Frameshift Position'].append("No Frameshift")
                    results['Segment 1'].append(sequence2)
                    results['Segment 2'].append(" ")
                    results['Gap'].append(gap_str)
                    results['Frameshift Direction'].append(
                        f'The sequence is present in {frame_name} without any frameshift'
                    )
                    results['Nucleotide Title'].append(fasta1_headers[b])
                    results['Nucleotide Sequence'].append(trimmed_seq)
                    results['Protein Title'].append(fasta2_headers[a])
                    results['Protein Sequence'].append(sequence2)
                    results['Truncation for Nucleotide Sequence'].append(str(truncated))

            # Create DataFrame for current nucleotide-peptide pair
            results_df = pd.DataFrame(results)
            results_all = pd.concat([results_all, results_df], ignore_index=True)

    return results_all


def process_sequences_helper(args):
    """Helper function for parallel processing"""
    (header1, sequence1, fasta2_headers, fasta2_sequences, max_gap, 
     speed_up_by_blast, blast_output, codon_table, 
     max_transcript_length, max_flanking_seq) = args
     
    return process_sequences(
        [header1], [sequence1], fasta2_headers, fasta2_sequences, max_gap,
        speed_up_by_blast, blast_output, codon_table, 
        max_transcript_length, max_flanking_seq
    )


# Replace the process_sequences_parallel function with this updated version
def process_sequences_parallel(fasta1_headers_all, fasta1_sequences_all, fasta2_headers,
                               fasta2_sequences, max_gap, speed_up_by_blast, blast_output,
                               n_jobs, codon_table, max_transcript_length, max_flanking_seq):
    """Parallel processing with human-readable progress information"""
    results_all = pd.DataFrame()
    args_list = [(header1, sequence1, fasta2_headers, fasta2_sequences, max_gap,
                  speed_up_by_blast, blast_output, codon_table,
                  max_transcript_length, max_flanking_seq)
                 for header1, sequence1 in zip(fasta1_headers_all, fasta1_sequences_all)]

    # Create custom progress bar
    pbar = tqdm(total=len(args_list), desc="Processing sequences", unit="seq")
    start_time = time.time()
    processed_count = 0
    processed_results = []

    try:
        # Process in parallel
        with multiprocessing.Pool(processes=n_jobs) as pool:
            for result in pool.imap_unordered(process_sequences_helper, args_list):
                processed_results.append(result)
                processed_count += 1

                # Calculate metrics
                elapsed = time.time() - start_time
                avg_speed = elapsed / processed_count
                total_minutes = elapsed / 60

                # Update progress bar with human-readable information
                pbar.set_postfix_str(
                    f"Total: {total_minutes:.1f} min for {processed_count} seqs, "
                    f"Speed: {avg_speed / 60:.1f} min/seq"
                )
                pbar.update(1)
    finally:
        pbar.close()
        total_time = time.time() - start_time

        # Print final summary
        print(f"  Sequences processed: {processed_count}")
        print(f"  Total time: {total_time / 60:.1f} minutes")
        print(f"  Average speed: {total_time / processed_count / 60:.1f} minutes per sequence")
        print(f"  ({processed_count / (total_time / 60):.2f} sequences per minute)")

    # Collect results
    for result in processed_results:
        results_all = pd.concat([results_all, result], ignore_index=True)

    return results_all

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = sequence[::-1]
    reverse_complement_seq = ''.join(complement.get(base, base) for base in reverse_seq)
    return reverse_complement_seq


def is_valid_fasta(file_path):
    try:
        with open(file_path, "r") as file:
            first_line = file.readline().strip()
            return first_line.startswith(">")
    except:
        return False


def clean_header(header):
    """Remove special characters from headers that might cause TSV issues"""
    return re.sub(r'[\t\n\r]', ' ', header)


def main():
    parser = argparse.ArgumentParser(
        description='ShiftSCAN - Identify chimeric sequences in nucleotide databases through frameshift analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--nucleotide_file', required=True,
                        help='Path to nucleotide FASTA file (genome/transcriptome/repeatome)')
    parser.add_argument('--peptide_file', required=True,
                        help='Path to query peptide FASTA file')

    # Processing parameters
    parser.add_argument('--max_gap', type=int, default=2,
                        help='Maximum allowed gap size in nucleotides (frameshift value)')
    parser.add_argument('--codon_table', type=int, default=1,
                        help='ID of the NCBI codon table (1=Standard Code)')
    parser.add_argument('--num_threads', type=int, default=1,
                        help='Number of parallel processing threads')

    # BLAST parameters
    parser.add_argument('--blast', action='store_true',
                        help='Enable BLAST pre-filtering for acceleration (not recommended)')
    parser.add_argument('--threshold', type=float, default=5.0,
                        help='BLAST word inclusion threshold (higher = slower)')
    parser.add_argument('--word_size', type=int, default=2,
                        help='BLAST word size parameter')
    parser.add_argument('--evalue', type=float, default=1000.0,
                        help='BLAST e-value threshold')

    # Output control
    parser.add_argument('--output', default='results',
                        help='Base name for the TSV output file')
    parser.add_argument('--max_transcript_length', type=int, default=30000,
                        help='Max length of a transcript sequence (nucleotides)')
    parser.add_argument('--max_flanking_seq', type=int, default=500,
                        help='Max length of flanking sequences to include (nucleotides)')

    # Analysis options
    parser.add_argument('--no_reverse_complement_check', action='store_true',
                        help='Disable reverse complement analysis')

    args = parser.parse_args()

    # Print thread information
    print(f"\nNumber of available threads: {args.num_threads}\n")

    # Set up codon table
    if args.codon_table == 0:
        codon_table = codon_table_default
        print("Using default codon table")
    else:
        try:
            codon_table = biotite.sequence.CodonTable.load(args.codon_table).codon_dict()
            print(f"Using codon table #{args.codon_table}")
        except Exception as e:
            sys.exit(f"Error loading codon table: {str(e)}")

    # Validate input files
    for path in [args.nucleotide_file, args.peptide_file]:
        if not os.path.exists(path):
            sys.exit(f"Error: File not found - {path}")
        if not is_valid_fasta(path):
            sys.exit(f"Error: Invalid FASTA format - {path}")

    # Validate output file name
    if any(c in args.output for c in r'\/:*?"<>|'):
        sys.exit("Error: Invalid characters in output name")

    # Read input files
    try:
        nuc_headers = []
        nuc_seqs = []
        pep_headers = []
        pep_seqs = []

        # Parse and clean headers
        raw_nuc_headers, raw_nuc_seqs = parse_fasta(args.nucleotide_file)
        for header, seq in zip(raw_nuc_headers, raw_nuc_seqs):
            nuc_headers.append(clean_header(header))
            nuc_seqs.append(seq)

        raw_pep_headers, raw_pep_seqs = parse_fasta(args.peptide_file)
        for header, seq in zip(raw_pep_headers, raw_pep_seqs):
            pep_headers.append(clean_header(header))
            pep_seqs.append(seq)

        print(f"Loaded {len(nuc_headers)} nucleotide sequences")
        print(f"Loaded {len(pep_headers)} peptide sequences")
    except Exception as e:
        sys.exit(f"Error reading FASTA files: {str(e)}")

    # BLAST processing (optional)
    blast_df = pd.DataFrame()
    if args.blast:
        print("\n🚀 Initializing BLAST acceleration...")
        try:
            blast_cmd = [
                'tblastn',
                '-query', args.peptide_file,
                '-subject', args.nucleotide_file,
                '-seg', 'no',
                '-soft_masking', 'false',
                '-word_size', str(args.word_size),
                '-threshold', str(args.threshold),
                '-evalue', str(args.evalue),
                '-outfmt', '6 qseqid sseqid pident',
                '-max_target_seqs', str(len(nuc_headers)),
                '-num_threads', str(args.num_threads)
            ]

            # Run BLAST and capture output
            result = subprocess.run(
                blast_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )

            # Process BLAST output
            blast_lines = result.stdout.strip().split('\n')
            if blast_lines and blast_lines[0]:
                blast_data = [line.split('\t') for line in blast_lines]
                blast_df = pd.DataFrame(blast_data)
                # Clean BLAST headers
                blast_df[0] = blast_df[0].apply(clean_header)
                blast_df[1] = blast_df[1].apply(clean_header)
                blast_df = blast_df[blast_df[2] == '100']  # Filter 100% identity matches
                print(f"✅ BLAST found {len(blast_df)} preliminary matches")
            else:
                print("⚠️ BLAST found no matches")

        except subprocess.CalledProcessError as e:
            sys.exit(f"BLAST failed: {e.stderr}")
        except FileNotFoundError:
            sys.exit("Error: tblastn not found. Please install BLAST+ and ensure it's in your PATH")

    # Main processing
    print("\n🔍 Analyzing sequences...")
    start_time = time.time()

    try:
        # Forward processing
        print("\nProcessing forward sequences:")
        forward_results = process_sequences_parallel(
            nuc_headers, nuc_seqs,
            pep_headers, pep_seqs,
            args.max_gap,
            args.blast,
            blast_df,
            args.num_threads,
            codon_table,
            args.max_transcript_length,
            args.max_flanking_seq
        )
        if not forward_results.empty:
            forward_results['Frame Direction'] = 'Forward'

        # Reverse complement processing
        reverse_results = pd.DataFrame()
        if not args.no_reverse_complement_check:
            print("\n🔄 Processing reverse complements:")
            rev_seqs = [reverse_complement(s) for s in nuc_seqs]
            reverse_results = process_sequences_parallel(
                nuc_headers, rev_seqs,
                pep_headers, pep_seqs,
                args.max_gap,
                args.blast,
                blast_df,
                args.num_threads,
                codon_table,
                args.max_transcript_length,
                args.max_flanking_seq
            )
            if not reverse_results.empty:
                reverse_results['Frame Direction'] = 'Reverse'

        # Combine results
        combined = pd.concat([forward_results, reverse_results], ignore_index=True)
        if combined.empty:
            print("⚠️ No results found")
        else:
            # Reorder columns to match specification
            column_order = [
                'Type', 'Frameshift Position', 'Segment 1', 'Segment 2', 'Gap',
                'Frameshift Direction', 'Nucleotide Title', 'Nucleotide Sequence',
                'Protein Title', 'Protein Sequence', 'Frame Direction',
                'Truncation for Nucleotide Sequence'
            ]
            combined = combined[column_order]

    except Exception as e:
        sys.exit(f"Processing failed: {str(e)}")

    # Generate TSV output
    output_tsv = f"{args.output}.tsv"
    try:
        if not combined.empty:
            # Save as TSV with proper formatting
            combined.to_csv(output_tsv, sep='\t', index=False)
            print(f"\n💾 Results saved to {output_tsv}")
            print(f"   Total entries: {len(combined)}")
        else:
            print("⚠️ No results to save")

    except Exception as e:
        sys.exit(f"Failed to write output: {str(e)}")

    # Save parameters
    params = f"""Analysis Parameters
{'=' * 18}
Nucleotide file: {args.nucleotide_file}
Peptide file: {args.peptide_file}
Codon table: {'Default' if args.codon_table == 0 else args.codon_table}
Max gap size: {args.max_gap}
Max transcript length: {args.max_transcript_length}
Max flanking sequence: {args.max_flanking_seq}
BLAST acceleration: {'Enabled' if args.blast else 'Disabled'}
BLAST threshold: {args.threshold}
BLAST word size: {args.word_size}
BLAST e-value: {args.evalue}
Reverse complement check: {'Disabled' if args.no_reverse_complement_check else 'Enabled'}
Number of threads: {args.num_threads}
Total sequences processed: {len(nuc_headers)} nucleotides, {len(pep_headers)} peptides
Processing time: {time.time() - start_time:.2f} seconds
"""
    try:
        with open(f"{args.output}_parameters.txt", 'w') as f:
            f.write(params)
        print(f"📝 Parameters saved to {args.output}_parameters.txt")
    except Exception as e:
        print(f"⚠️ Failed to save parameters: {str(e)}")

def run():
    main()

if __name__ == "__main__":
    run()
