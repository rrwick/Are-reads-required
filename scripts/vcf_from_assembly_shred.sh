#!/usr/bin/env bash

# This script produces a VCF from an assembly file by shredding the assembly into fake reads and
# then calling variants with a read-based approach (BWA and Freebayes).

# Note that the number of read pairs (962200) is set to give 100x depth for the S. aureus NRS384
# genome, and it may not be ideal for genomes of different sizes.

# Requirements:
# - wgsim
# - bwa
# - samtools
# - freebayes
# - bcftools 


# Get arguments.
ref=$1
assembly=$2
threads=$3

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Copy the reference genome to the temporary directory and build indices.
cp "$ref" "$temp_dir"/ref.fasta
bwa index "$temp_dir"/ref.fasta
samtools faidx "$temp_dir"/ref.fasta

# Simulate error-free reads at 100x depth.
wgsim -e 0 -r 0 -1 150 -2 150 -N 962200 "$assembly" "$temp_dir"/1.fastq "$temp_dir"/2.fastq

# Align the fake reads to the reference.
bwa mem -t "$threads" "$temp_dir"/ref.fasta "$temp_dir"/1.fastq "$temp_dir"/2.fastq | samtools sort > "$temp_dir"/alignments.bam

# Call SNPs with freebayes
freebayes -f "$temp_dir"/ref.fasta --haplotype-length -1 -m 10 -q 10 -p 1 --min-coverage 2 "$temp_dir"/alignments.bam | bcftools view -i 'QUAL>=100'
