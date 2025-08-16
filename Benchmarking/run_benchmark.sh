#!/bin/bash
#SBATCH --time=2-0:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --exclusive
 
set -euo pipefail
if [[ "${DEBUG:-0}" == "1" ]]; then set -x; fi

# --- Config ---
cd benchmarking/

ROOT_DIR="$(pwd)"
NUC_DIR="cDNA_sequences"
PEP_DIR="chimeric_peptides"
OUT_ROOT="bench_out"
LOG_CSV="benchmark_results.csv"
THREADS=$(seq 1 128)

ABS_NUC_DIR="$ROOT_DIR/$NUC_DIR"
ABS_PEP_DIR="$ROOT_DIR/$PEP_DIR"

mkdir -p "$OUT_ROOT"

# --- Helpers ---
species_from_base () { printf "%s" "${1%%.*}"; }

# Replace dots with underscores, then allow only [A-Za-z0-9_-]
sanitize_name () {
  local s="$1"
  s="${s//./_}"
  s="$(printf "%s" "$s" | sed 's/[^A-Za-z0-9_-]/_/g')"
  s="$(printf "%s" "$s" | sed 's/___*/_/g')"
  s="$(printf "%s" "$s" | sed 's/^_//; s/_$//')"
  printf "%s" "$s"
}

trap 'rc=$?; echo "ERROR: Script aborted (exit $rc)"; exit $rc' ERR

# CSV header
if [[ ! -f "$LOG_CSV" ]]; then
  echo "species,threads,elapsed_seconds,start_iso,end_iso,output_dir,command,nucleotide_file,peptide_file" > "$LOG_CSV"
fi

# Info about available CPUs (best-effort)
max_avail="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}"
echo "Note: this job/node seems to have up to $max_avail cores available."

# Main loop over full nucleotide files
for nucl_rel in "$NUC_DIR"/*cdna.all.fa; do
  base="$(basename "$nucl_rel" .fa)"   # e.g., Homo_sapiens.GRCh38.cdna.all
  pep_rel="$PEP_DIR/chimeric_peptides_${base}_combined.fasta"

  # Absolute input paths (so we can cd elsewhere safely)
  nucl="$ABS_NUC_DIR/$(basename "$nucl_rel")"
  pep="$ABS_PEP_DIR/$(basename "$pep_rel")"

  if [[ ! -s "$pep" ]]; then
    echo "Skipping $base: missing peptide file -> $pep" >&2
    continue
  fi
  if [[ ! -s "$nucl" ]]; then
    echo "Skipping $base: missing nucleotide file -> $nucl" >&2
    continue
  fi

  species="$(species_from_base "$base")"
  safe_base="$(sanitize_name "$base")"

  echo "== Pair =="
  echo "  nucl: $nucl"
  echo "  pep : $pep"
  echo "  base: $base"
  echo "  species: $species"
  echo "  safe_base: $safe_base"

  for t in $THREADS; do
    out_dir="$OUT_ROOT/${species}/full/t${t}"
    mkdir -p "$out_dir"

    if (( t > max_avail )); then
      echo "  [warn] threads=$t > available=$max_avail (run will still proceed)"
    fi

    # Build command using ABSOLUTE input paths; run from out_dir so --output has no slashes
    cmd=(shiftscan
         --nucleotide_file "$nucl"
         --peptide_file    "$pep"
         --output          "$safe_base"
         --num_threads     "$t"
         --no_reverse_complement_check) 

    echo ">> ${species} | threads=${t}"
    echo "   out_dir: $out_dir"
    echo "   command: (cd \"$out_dir\" && ${cmd[*]})"

    start_epoch=$(date +%s)
    start_iso=$(date -Is)

    ( cd "$out_dir" && "${cmd[@]}" )

    end_epoch=$(date +%s)
    end_iso=$(date -Is)
    elapsed=$(( end_epoch - start_epoch ))

    printf "%s,%d,%d,%s,%s,%s,%q,%s,%s\n" \
      "$species" "$t" "$elapsed" "$start_iso" "$end_iso" "$out_dir" "${cmd[*]}" "$nucl" "$pep" \
      >> "$LOG_CSV"
  done
done

echo "Done. Results: $LOG_CSV"
echo "Outputs: $OUT_ROOT/<species>/full/t<threads>/"
