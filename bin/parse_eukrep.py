#!/usr/bin/env python

"""
Parse EukRep

Last updated 7/11/22
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import gzip
import argparse
import pandas as pd

from collections import defaultdict
from Bio import SeqIO

def euk_parser(ori_fasta, eukrep_out, stb, s2l=None):
    # Get scaffold2length
    if s2l is None:
        s2l = get_s2l(ori_fasta)
    
    # Figure out which scaffolds are eukaryotic
    euk_scaffs = get_scaff_list(eukrep_out)
    
    # Load the stb
    stb = load_stb(stb)
    
    # Calculate % euk
    edb = calc_p_euk(euk_scaffs, stb, s2l)
    
    return edb
    
def get_s2l(ori_fasta):
    # Get genome to length
    if ori_fasta.endswith('.gz'):
        handle = gzip.open(ori_fasta, "rt")
    else:
        handle = open(ori_fasta, "r")
        
    scaff2sequence = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    s2l = {s: len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}
    
    handle.close()
    
    return s2l

def get_scaff_list(eukrep_out):
    if eukrep_out.endswith('.gz'):
        handle = gzip.open(eukrep_out, "rt")
    else:
        handle = open(eukrep_out, "r")
        
    scaffs = []
    for line in handle.readlines():
        if line.startswith('>'):
            scaffs.append(line[1:].split()[0].strip())
            
    handle.close()
    return scaffs

def load_stb(stb_file):
    if stb_file.endswith('.gz'):
        handle = gzip.open(stb_file, "rt")
    else:
        handle = open(stb_file, "r")
        
    stb = {}
    for line in handle.readlines():
        stb[line.split()[0].strip()] = line.split()[1].strip()
    
    handle.close()
    return stb

def calc_p_euk(euk_scaffs, stb, s2l):
    b2s = defaultdict(list)
    for s, b in stb.items():
        b2s[b].append(s)
    
    table = defaultdict(list)
    euk_bins = set([stb[s] for s in euk_scaffs if s in stb])
    for b in euk_bins:
        scaffs = set(b2s[b])
        e_scaffs = set(euk_scaffs).intersection(scaffs)
        
        tlen = sum([s2l[s] for s in scaffs])
        elen = sum([s2l[s] for s in e_scaffs])
        
        table['bin'].append(b)
        table['scaffolds'].append(len(scaffs))
        table['length'].append(tlen)
        table['euk_scaffolds'].append(len(e_scaffs))
        table['euk_length'].append(elen)
    
    edb = pd.DataFrame(table)
    if len(edb) > 0:
        edb['p_euk_scaffs'] = (edb['euk_scaffolds'] / edb['scaffolds']) * 100
        edb['p_euk_length'] = (edb['euk_length'] / edb['length']) * 100
        edb = edb.sort_values('p_euk_length', ascending=False).reset_index(drop=True)
    
    return edb
                

def main(args):
    print(args)
    edb = euk_parser(args.fasta, args.eukrep, args.stb)
    print(edb)
    edb.to_csv(args.out, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta", action="store", help='Location of .fasta file EukRep was run on')
    parser.add_argument("-e", "--eukrep", action="store", help='Location of EukRep output')
    parser.add_argument("-s", "--stb", action="store", help='Location of stb file for parsing')
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
