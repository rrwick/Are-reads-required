#!/usr/bin/env python3
"""
This script takes in a FASTA MSA and it drops all sites that are identical. Each character is
treated as unique, regardless of whether it is a DNA base or not, so gaps will be included in the
final output. The reduced FASTA MSA will be outputted to stdout.

It also counts the number of total pairwise differences (all possible pairwise combinations), and
outputs this to stderr. Two counts are given: the first ignores gaps and the second counts gaps.
"""

import gzip
import itertools
import pytest
import sys


def main():
    headers, seqs = load_input_alignment(sys.argv[1])
    new_seqs = drop_invariant_sites(seqs)
    print_alignment(headers, new_seqs)
    total_no_gaps, total_with_gaps = total_pairwise_difference_count(new_seqs)
    print(f'{total_no_gaps}\t{total_with_gaps}', file=sys.stderr)


def load_input_alignment(filename):
    headers, seqs = zip(*iterate_fasta(filename))
    seq_len = len(seqs[0])
    assert all(len(seq) == seq_len for seq in seqs), "All sequences must be the same length"
    return list(headers), list(seqs)


def drop_invariant_sites(seqs):
    seq_len = len(seqs[0])
    new_seqs = [[] for _ in seqs]
    for i in range(seq_len):
        if all(s[i] == seqs[0][i] for s in seqs):  # all seqs are identical at position i
            continue
        for seq, new_seq in zip(seqs, new_seqs):
            new_seq.append(seq[i])
    return [''.join(s) for s in new_seqs]


def print_alignment(headers, seqs):
    for header, seq in zip(headers, seqs):
        print(header)
        print(seq)


def total_pairwise_difference_count(seqs):
    total_no_gaps, total_with_gaps = 0, 0
    for seq_a, seq_b in itertools.combinations(seqs, 2):
        total_no_gaps += difference_count(seq_a, seq_b, count_gaps=False)
        total_with_gaps += difference_count(seq_a, seq_b, count_gaps=True)
    return total_no_gaps, total_with_gaps


def difference_count(seq_a, seq_b, count_gaps):
    assert len(seq_a) == len(seq_b)
    return sum(a != b and (count_gaps or (a != '-' and b != '-'))
               for a, b in zip(seq_a, seq_b))


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    return gzip.open if get_compression_type(filename) == 'gz' else open


def iterate_fasta(filename):
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        header, sequence = '', []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if header:
                    yield header, ''.join(sequence)
                header, sequence = line[1:], []
            else:
                sequence.append(line.upper())
        if header:
            yield header, ''.join(sequence)


if __name__ == '__main__':
    main()


# Automated tests with Pytest

def test_difference_count():
    assert difference_count('', '', False) == 0
    assert difference_count('', '', True) == 0

    assert difference_count('ACGT', 'ACGT', False) == 0
    assert difference_count('ACGT', 'ACGT', True) == 0
    
    assert difference_count('ACGT', 'AC-T', False) == 0
    assert difference_count('ACGT', 'AC-T', True) == 1
    
    assert difference_count('ACGT', 'TGCA', False) == 4
    assert difference_count('ACGT', 'TGCA', True) == 4
    
    assert difference_count('ACGT----', '----ACGT', False) == 0
    assert difference_count('ACGT----', '----ACGT', True) == 8


def test_total_pairwise_difference_count():
    assert total_pairwise_difference_count(['', '', '', '']) == (0, 0)

    assert total_pairwise_difference_count(['ACGT', 'ACGT', 'ACGT', 'ACGT']) == (0, 0)

    assert total_pairwise_difference_count(['AAGT', 'ACGT', 'ACGT', 'ACGT']) == (3, 3)
    assert total_pairwise_difference_count(['ACGT', 'ACCT', 'ACGT', 'ACGT']) == (3, 3)
    assert total_pairwise_difference_count(['ACGT', 'ACGT', 'ACGG', 'ACGT']) == (3, 3)
    assert total_pairwise_difference_count(['ACGT', 'ACGT', 'ACGT', 'TCGT']) == (3, 3)

    assert total_pairwise_difference_count(['-CGT', 'ACGT', 'ACGT', 'ACGT']) == (0, 3)
    assert total_pairwise_difference_count(['ACGT', 'A-GT', 'ACGT', 'ACGT']) == (0, 3)
    assert total_pairwise_difference_count(['ACGT', 'ACGT', 'AC-T', 'ACGT']) == (0, 3)
    assert total_pairwise_difference_count(['ACGT', 'ACGT', 'ACGT', 'ACG-']) == (0, 3)

    assert total_pairwise_difference_count(['ACGA', 'AGG-', 'AC-T', 'ATGT']) == (7, 13)
