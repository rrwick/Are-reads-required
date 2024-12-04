#!/usr/bin/env bash

# This script produces a VCF from an assembly file using SKA.

# In order to produce a VCF compatible with vcfdist, it adds reference sequence lenghts to the
# header, so as written, this script is only compatible with my specific reference genome
# (S. aureus NRS384).

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

ska build -o "$temp_dir"/ska.skf -k 31 "$assembly"
ska weed --filter no-ambig "$temp_dir"/ska.skf
ska map "$ref" "$temp_dir"/ska.skf -f vcf | bcftools view -e'ALT="."' | sed 's/##contig=<ID=chromosome>/##contig=<ID=chromosome,length=2879034>/' | sed 's/##contig=<ID=plasmid_1>/##contig=<ID=plasmid_1,length=4439>/' | sed 's/##contig=<ID=plasmid_2>/##contig=<ID=plasmid_2,length=3125>/'
