#!/usr/bin/env bash

# This script produces a VCF from an assembly file using SKA.

# Requirements:
# - ska (the newer Rust-based version at github.com/bacpop/ska.rust)
# - bcftools

# Get arguments.
ref=$1
assembly=$2

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

samtools faidx "$ref"
ska build -o "$temp_dir"/ska.skf -k 31 "$assembly"
ska weed --filter no-ambig "$temp_dir"/ska.skf
ska map "$ref" "$temp_dir"/ska.skf -f vcf | bcftools view -e'ALT="."' | bcftools reheader -f "${ref}.fai"
