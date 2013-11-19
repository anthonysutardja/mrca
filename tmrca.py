#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""


def read_fasta_sequence_to_str(filename):
    """
    Reads a FASTA formatted file to a string.

    Returns a tuple containing the sequence's name and the sequence.
    """
    with open(filename) as f:
        lines = [line.strip() for line in f.readlines()]
        sequence_name = lines.pop(0).replace('>', '')
        return (sequence_name, ''.join(lines))


def main():
    pass


if __name__ == '__main__':
    main()
