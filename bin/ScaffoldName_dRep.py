#!/usr/bin/env python

# Matt Olm
# mattolm@stanford.edu
# 08.01.22
"""
Create a copy of a .fasta file where no scaffolds are repeated
"""
import os
import gzip
import textwrap
import argparse

from Bio import SeqIO

__author__ = "Matt Olm"
__license__ = "MIT"
__version__ = "0.1.0"

def main(args):
    name_drep(args.input, args.output)

def name_drep(input, output):
    IDs = set()
    Repeats = set()

    nw = open(output, 'a')

    if input[-3:] == '.gz':
        handle = gzip.open(input, "rt")
        seqs = SeqIO.parse(handle, "fasta")

    else:
        seqs = SeqIO.parse(input, "fasta")

    for seq_record in seqs:
        id = str(seq_record.id).strip()
        if id in IDs:
            Repeats.add(id)
        else:
            nw.write('\n'.join([">{0}".format(id), str(seq_record.seq), '']))
        IDs.add(id)

    nw.close()
    if input[-3:] == '.gz':
        handle.close()

    print(f"{len(Repeats)} scaffolds are repeated")

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         ScaffoldName_dRep creates a copy of a .fasta file where no scaffolds are repeated
         '''))

    parser.add_argument("-i", "--input", help='input .fasta file')
    parser.add_argument("-o", "--output", help='output .fasta file')

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
