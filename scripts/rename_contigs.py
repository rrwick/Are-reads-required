#!/usr/bin/env python3
"""
This script takes in a polished completed assembly (FASTA format) and renames the contigs,
overwriting the file. It makes the following assumptions:
* The FASTA file has one line per sequence (no multi-line sequences).
* The first contig is the chromosome and the rest are plasmids.
* All contigs are circular.
"""

import sys


def main():
    plasmid_num = 0
    output_lines = []
    with open(sys.argv[1], 'rt') as f:
        for header in f:
            seq = next(f).rstrip('\n')
            if plasmid_num == 0:
                name = 'chromosome'
            else:
                name = f'plasmid_{plasmid_num}'
            plasmid_num += 1
            output_lines.append(f'>{name} circular=true length={len(seq)}')
            output_lines.append(seq)
    with open(sys.argv[1], 'wt') as f:
        for line in output_lines:
            f.write(line)
            f.write('\n')


if __name__ == '__main__':
    main()
