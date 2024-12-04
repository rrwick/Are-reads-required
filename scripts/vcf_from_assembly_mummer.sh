#!/usr/bin/env bash

# This script produces a VCF from an assembly file using MUMmer.

# Since the all2vcf output isn't quite compatible with vcfdist, I manipulate the result a bit by
# adding a 'reference' sample column.

# Requirements:
# - MUMmer v4 (specifically the dnadiff script)
# - all2vcf (github.com/MatteoSchiavinato/all2vcf)

# Get arguments.
ref=$1
assembly=$2

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

dnadiff -p "$temp_dir"/mummer "$ref" "$assembly"
all2vcf mummer --snps "$temp_dir"/mummer.snps --reference "$ref" --output-header --no-Ns | awk 'BEGIN{OFS="\t"} /^##/ {print $0} /^#CHROM/ {print $0, "FORMAT", "unknown"} !/^#/ {NF--; print $0, "GT", "1"}'
