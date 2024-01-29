#!/usr/bin/env python

"""
Parse coverM

Last updated 8/24/23
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import sys
import site
import argparse
import importlib
import subprocess

from collections import defaultdict

def process_coverm(db, s2b, bin2length):
    """
    s2b = scaffold 2 bin
    bin2length = bin 2 length
    """
    # Add bin
    ndb = db.copy()
    ndb['bin'] = ndb['Contig'].map(s2b)
    #assert len(ndb[ndb['bin'].isna()]) == 0, ndb[ndb['bin'].isna()]

    # Loop through bins
    table = defaultdict(list)
    assert len(ndb['Sample'].unique()) == 1
    sample = ndb['Sample'].iloc[0]
    for bin, bdb in ndb.groupby('bin'):
        true_length = bin2length[bin]
        
        covered_length = bdb['Covered Bases'].sum()

        bdb['ml'] = bdb['Mean'] * bdb['Length']
        coverage = bdb['ml'].sum() / true_length

        table['sample'].append(sample)
        table['genome'].append(bin)
        table['breadth'].append(covered_length/true_length)
        table['coverage'].append(coverage)
        table['mapped_reads'].append(bdb['Read Count'].sum())
        table['length'].append(true_length)
    
    return table

def generate_bin_length_dict(scaffold_bin_dict, scaffold_length_dict):
    bin_length_dict = {}

    for scaffold, bin_id in scaffold_bin_dict.items():
        length = scaffold_length_dict.get(scaffold)
        if length is not None:
            bin_length_dict[bin_id] = bin_length_dict.get(bin_id, 0) + length

    return bin_length_dict

def generate_scaffold_length_dict(fasta_file):
    scaffold_length_dict = {}
    current_scaffold = None
    current_length = 0

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('>'):
                if current_scaffold is not None:
                    scaffold_length_dict[current_scaffold] = current_length
                    current_length = 0

                current_scaffold = line[1:].split()[0]
            else:
                current_length += len(line)

        # Add the length of the last scaffold
        if current_scaffold is not None:
            scaffold_length_dict[current_scaffold] = current_length

    return scaffold_length_dict

def install_and_import(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", package])

    # Add user site-packages to sys.path
    user_site = site.getusersitepackages()
    if user_site not in sys.path:
        sys.path.append(user_site)
    importlib.reload(site)
    __import__(package)

def main(args):
    install_and_import('pandas')
    import pandas as pd

    odb = pd.read_csv(args.coverm, sep='\t')
    
    if len(args.stb) == 0:
        print('just returning original')
        db = odb

    else:
        stb = pd.read_csv(args.stb[0], sep='\t', names=['scaffold', 'bin']).set_index('scaffold')['bin'].to_dict()
        s2l = generate_scaffold_length_dict(args.fasta)
        bin2length = generate_bin_length_dict(stb, s2l)

        table = process_coverm(odb, stb, bin2length)
        db = pd.DataFrame(table)

    db.to_csv(args.out, index=False, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--coverm", action="store", help='Location of coverm .tsv output')
    parser.add_argument("-f", "--fasta", action="store", help='Location of fasta file')
    parser.add_argument("-s", "--stb", action="store", help='Location of stb file', required=False, nargs='*')
    parser.add_argument("-o", "--out", action="store", help='Name of the output table to create')

    # # Optional argument flag which defaults to False
    # parser.add_argument("-f", "--flag", action="store_true", default=False)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
