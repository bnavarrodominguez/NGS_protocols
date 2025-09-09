#!/usr/bin/env bash
set -euo pipefail

# Usage: fastp_like_trimmomatic.sh <LIB> <suffix: fastq.gz|fq.gz> <threads> <adapters.fa|auto>
# Example: ./fastp_like_trimmomatic.sh SAMPLE fastq.gz 8 adapters.fa
#          ./fastp_like_trimmomatic.sh SAMPLE fq.gz     8 auto

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <LIB> <suffix> <threads> <adapters.fa|auto>" >&2
  exit 1
fi

lib="$1"
suf="$2"
thr="$3"
adapters="$4"

r1="${lib}_1.${suf}"
r2="${lib}_2.${suf}"

# deps
command -v fastp  >/dev/null 2>&1 || { echo "Error: fastp not found in PATH"  >&2; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo "Error: fastqc not found in PATH" >&2; exit 1; }

[[ -s "$r1" && -s "$r2" ]] || { echo "Error: missing inputs: $r1 / $r2" >&2; exit 1; }

out1="${lib}_paired_1.fastq.gz"
out2="${lib}_paired_2.fastq.gz"
un1="${lib}_unpaired_1.fastq.gz"
un2="${lib}_unpaired_2.fastq.gz"
json="${lib}.fastp.json"
html="${lib}.fastp.html"

# Build adapter-trimming arguments
adapter_args=()
if [[ "$adapters" == "auto" ]]; then
  adapter_args+=( --detect_adapter_for_pe )
else
  [[ -s "$adapters" ]] || { echo "Error: adapters FASTA not found: $adapters" >&2; exit 1; }
  adapter_args+=( --adapter_fasta "$adapters" )
fi

echo "[fastp] Trimming (adapter + sliding window Q20 + MINLEN 20)..."
fastp \
  -i "$r1" -I "$r2" \
  -o "$out1" -O "$out2" \
  "${adapter_args[@]}" \
  --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 20 \
  --unpaired1 "$un1" --unpaired2 "$un2" \
  --trim_poly_g \
  -j "$json" -h "$html" -w "$thr"

echo "[fastqc] Checking trimmed outputs..."
fastqc -t "$thr" "$out1" "$out2"

echo "Done."
echo "Outputs:"
echo "  Paired    : $out1  $out2"
echo "  Unpaired  : $un1   $un2"
echo "  Reports   : $html  ($json)"

