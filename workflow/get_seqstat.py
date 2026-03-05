#!/usr/bin/env python3

import sys
import os
from statistics import median

def read_fasta_lengths(fasta_file):
    lengths = []
    seq_len = 0

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line)

        if seq_len > 0:
            lengths.append(seq_len)

    return lengths


def process_file(fasta):
    lengths = read_fasta_lengths(fasta)

    if lengths:
        n = len(lengths)
        med = int(median(lengths))
    else:
        n = 0
        med = 0

    basename = os.path.splitext(os.path.basename(fasta))[0]
    print(f"{basename}\t{n}\t{med}")


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: get_input.py <fasta1> [fasta2 ...]\n")
        sys.exit(1)

    for fasta in sys.argv[1:]:
        process_file(fasta)


if __name__ == "__main__":
    main()