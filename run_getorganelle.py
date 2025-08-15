#!/usr/bin/env python3

import sys
import subprocess
import os
import shutil
from glob import glob
from Bio import SeqIO

def usage():
    print("Usage: run_getorganelle.py NumberOfReads ListOfFiles Reference Prefix Threads ")
    sys.exit(1)

# Parse command-line arguments
if len(sys.argv) != 6:
    usage()

nreads = sys.argv[1]
list_file = sys.argv[2]
ref = sys.argv[3]
prefix = sys.argv[4]
thr = sys.argv[5]

output_dir = f"{prefix}_results"
os.makedirs(output_dir, exist_ok=True)

# Check list file exists
if not os.path.isfile(list_file):
    print(f"Error: File '{list_file}' not found.")
    sys.exit(1)

# Read file list (one pair per two lines)
with open(list_file) as f:
    pairs = [line.strip() for line in f if line.strip()]

if len(pairs) % 2 != 0:
    print("Error: List file should contain an even number of lines (R1 and R2 pairs).")
    sys.exit(1)

# Loop over pairs
for i in range(0, len(pairs), 2):
    r1 = pairs[i]
    r2 = pairs[i + 1]

    # Derive sample name
    basename = os.path.basename(r1).split('.')[0].rstrip('_R1').rstrip('_1')
    foldername = f"{prefix}_{basename}"

    print(f"\n>>> Running GetOrganelle for {basename}")

    cmd = [
        "get_organelle_from_reads.py",
        "-1", r1,
        "-2", r2,
        "-o", foldername,
        "-s", ref,
        "--genes", ref,
        "-R", "10",
        "-k", "21,45,65,85,105",
        "--max-reads", nreads,
        "-t", thr,
        "-F", "anonym"
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f">>> Finished {basename}")

        # Find all FASTA files in the output directory
        fasta_candidates = glob(os.path.join(foldername, "*.fasta"))
        gfa_candidates = glob(os.path.join(foldername, "*.gfa"))

        if fasta_candidates:
            dest_fasta = os.path.join(output_dir, f"{basename}.fasta")
            total = 0
            with open(dest_fasta, "w") as out_f:
                for fasta in fasta_candidates:
                    for record in SeqIO.parse(fasta, "fasta"):
                        record.id = f"{basename}_{record.id}"
                        record.description = ""
                        SeqIO.write(record, out_f, "fasta")
                        total += 1
            print(f"✅Wrote {total} sequences to {dest_fasta}")
        else:
            print(f" No FASTA files found for {basename}")

        if gfa_candidates:
            dest_gfa = os.path.join(output_dir, f"{basename}.gfa")
            shutil.copyfile(gfa_candidates[0], dest_gfa)
        else:
            print(f" No GFA file found for {basename}")

    except subprocess.CalledProcessError as e:
        print(f" Error running GetOrganelle for {basename}")
        print(e)

