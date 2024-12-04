#!/usr/bin/env python3
"""
This script takes in minimap2 alignments (query=assembly, target=reference) and outputs three
numbers:
1) how many bases in the reference were missed (not covered by any alignments)
2) how many bases in the reference were covered once
3) how many bases in the reference were covered multiple times
"""

import fileinput


def main():
    alignments = []
    for line in fileinput.input():
        alignments.append(Alignment(line))

    ref_lengths = {'chromosome': 2879034, 'plasmid_1': 4439, 'plasmid_2': 3125}

    alignments_by_ref = {r: [] for r in ref_lengths.keys()}
    for a in alignments:
        assert a.ref_length == ref_lengths[a.ref_name]
        alignments_by_ref[a.ref_name].append(a)

    missing_total, once_total, multiple_total = 0, 0, 0
    for ref_name, alignments in alignments_by_ref.items():
        ref_length = ref_lengths[ref_name]
        depths = [0] * ref_lengths[ref_name]
        for a in alignments:
            for i in range(a.ref_start, a.ref_end):
                depths[i] += 1
        missing = depths.count(0)
        once = depths.count(1)
        multiple = ref_length - missing - once
        missing_total += missing
        once_total += once
        multiple_total += multiple

    assert missing_total + once_total + multiple_total == sum(ref_lengths.values())
    print(f'{missing_total}\t{once_total}\t{multiple_total}')


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])


if __name__ == '__main__':
    main()
