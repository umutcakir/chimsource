import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

def open_fasta(fasta_file):
    """Open fasta or gzipped fasta file."""
    if fasta_file.endswith(".gz"):
        return gzip.open(fasta_file, "rt")
    else:
        return open(fasta_file, "r")

def get_chimeric_peptide(seq, mid_nt, shift):
    """
    seq: full nucleotide sequence (string)
    mid_nt: midpoint nucleotide index where frame shift occurs
    shift: frame shift integer (-2, -1, 0, +1, +2)
    returns: 10 amino acid peptide string (5 left + 5 right)
    """
    aa_each_side = 6
    nt_each_side = aa_each_side * 3

    # Left side (no frame shift)
    left_start = max(mid_nt - nt_each_side, 0)
    left_seq = seq[left_start:mid_nt]
    left_aa = str(Seq(left_seq).translate(to_stop=False))

    # Right side (frame shifted)
    right_start = mid_nt + shift
    right_seq = seq[right_start:right_start + nt_each_side]
    right_aa = str(Seq(right_seq).translate(to_stop=False))

    peptide = left_aa + right_aa

    # Validate length and absence of stop codons
    if len(peptide) < aa_each_side * 2:
        return None
    if "*" in peptide:
        return None

    return peptide[:aa_each_side * 2]

def translate_chimeric_peptides(fasta_path, output_prefix="chimeric_peptides", max_transcripts=200):
    with open_fasta(fasta_path) as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))

    if not records:
        print(f"No records found in: {fasta_path}")
        return

    output_dir = "simulated_chimeric_peptides"
    os.makedirs(output_dir, exist_ok=True)

    filename_with_ext = os.path.basename(fasta_path)
    if filename_with_ext.endswith(".fa.gz"):
        output_name = filename_with_ext[:-6]
    elif filename_with_ext.endswith(".fasta.gz"):
        output_name = filename_with_ext[:-9]
    else:
        output_name = os.path.splitext(filename_with_ext)[0]

    shifts = {
        "0_frame": 0,
        "+1_frame": 1,
        "+2_frame": 2,
        "-1_frame": -1,
        "-2_frame": -2
    }

    combined_output_file = os.path.join(
        output_dir, f"{output_prefix}_{output_name}_combined.fasta"
    )

    processed_count = 0

    with open(combined_output_file, "w") as out_f:
        for seq_record in records:
            if processed_count >= max_transcripts:
                break

            full_seq = str(seq_record.seq)
            mid_nt = len(full_seq) // 2

            peptides_for_transcript = []
            all_valid = True
            for label, shift in shifts.items():
                peptide = get_chimeric_peptide(full_seq, mid_nt, shift)
                if peptide:
                    peptides_for_transcript.append((label, peptide))
                else:
                    all_valid = False
                    break  # Skip transcript if any shift fails

            if all_valid and len(peptides_for_transcript) == 5:
                for label, peptide in peptides_for_transcript:
                    out_f.write(f">{seq_record.id}_{label}\n{peptide}\n")
                processed_count += 1

    print(
        f"Finished {output_name}: {processed_count} transcripts processed "
        f"({processed_count*5} peptides written) to {combined_output_file}"
    )

# Run for all species fasta files in input directory
input_dir = "cDNA_sequences"
for fasta_file in os.listdir(input_dir):
    fasta_path = os.path.join(input_dir, fasta_file)
    translate_chimeric_peptides(fasta_path, max_transcripts=20)
