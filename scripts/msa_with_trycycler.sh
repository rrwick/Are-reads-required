#!/usr/bin/env bash

# This script uses Trycycler to produce a whole genome MSA, which it outputs to stdout.

# Get arguments.
unaligned_fasta=$1
threads=$2

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Create a fake Trycycler cluster and run Trycycler MSA.
cp "$unaligned_fasta" "$temp_dir"/2_all_seqs.fasta
trycycler msa -c "$temp_dir" -t "$threads"

# Output the aligned FASTA.
cat "$temp_dir"/3_msa.fasta
