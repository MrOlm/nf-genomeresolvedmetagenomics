#!/usr/bin/env python

"""
Prepare dRep

Last updated 7/13/22
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

from dataclasses import dataclass
import os
import gzip
import glob
import shutil
import argparse
import pandas as pd

from collections import defaultdict

def main(args):
    print(args)

    # Figure out how many sets each folder has
    folders = [str(f) + '/' if f[-1] != '/' else f for f in list(args.folders)]
    sets = []
    for f in folders:
        for d in glob.glob(f + '*/'):
            sets.append(os.path.basename(d[:-1]))
    sets = set(sets)

    # Make a merged dRep folder for each set
    for s in sets:
        sname = f"{args.outbase}.{s}"

        dbs = []
        gf = f"{sname}_bins/"
        set_genomes = []
        database_genomes = []
        
        for f in folders:
            database_folder = 'database' in f

            # Make the genome info table
            gl = f"{f}{s}/{s}_genomeInfo.csv"
            if not os.path.isfile(gl):
                continue
            db = pd.read_csv(gl)
            dbs.append(db)

            # Make the folder of bins
            if not os.path.isdir(gf): os.mkdir(gf)
            genomes = glob.glob(f + s + '/bins/*')
            for g in genomes:
                set_genomes.append(g)
                shutil.copy2(g, gf + os.path.basename(g))
                if database_folder:
                    database_genomes.append(os.path.basename(g))
        database_genomes = set(database_genomes)


        # Save the genome info table
        gdb = pd.concat(dbs).reset_index(drop=True)
        gdb.to_csv(f"{sname}_genomeInfo.csv", index=False)

        # Make the extra scoring table
        assert len(set_genomes) == len(gdb), [len(set_genomes), len(gdb)]
        with open(f"{sname}_extraScore.tsv", 'w') as o:
            for g in set_genomes:
                bg = os.path.basename(g)
                if bg not in database_genomes:
                    o.write(f"{bg}\t{args.extra_score}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--folders"       , action="store",     help='Location of folders', nargs='*')
    parser.add_argument("-o", "--outbase"       , action="store",     help='Name of the output basename')

    parser.add_argument("-e", "--extra_score"   , action="store",     help='Extra score for genomes', default=0)

    # # Optional argument flag which defaults to False
    # parser.add_argument("-f", "--flag", action="store_true", default=False)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
