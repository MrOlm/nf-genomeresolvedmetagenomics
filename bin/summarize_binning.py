#!/usr/bin/env python

"""
Summarize Binning

v2.0 - 8/19/22
# In this version, phage / plasmid scaffolds in bacterial bins are NOT on their own in the .stb

v1.0 - 7/12/22
# NOTE: In this version, scaffolds that are phage / plasmid are removed from bacterial bins

(base) mattolm@mac-nugget:/LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin$ ~/user_data/Nextflow/nf-core-genomeresolvedmetagenomics/bin/summarize_binning.py -f /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Assemble/results_v1/assembly/MD_02_MEGAHIT.fasta.gz -s /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/parsestb/MD_02.stb -j /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/metabat2/MD_02.txt.gz -e /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/eukrep/MD_02.eukrep.csv -g /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/gtdbtk/gtdbtk.MD_02*summary.tsv -chm /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/checkm/MD_02.tsv -c /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v2/checkv/MD_02_results/complete_genomes.tsv -ch /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v2/checkv/MD_02_quality_summary.tsv -a /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v2/abricate/MD_02.txt -ts /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v2/trep/MD_02.trep.tax_fullScaffoldTaxonomy.tsv.gz -tg /LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/Wastewater/Nextflow/Bin/results_v3/trep/MD_02.trep.tax_fullGenomeTaxonomy.tsv.gz -o MD_02
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import os
import gzip
import argparse
import subprocess
import pandas as pd
import numpy as np

pd.set_option('mode.chained_assignment', None)

from collections import defaultdict
from Bio import SeqIO

def main(args):
    print(args)

    # 1) Create phage / plasmid table
    PPdb = phage_plasmid_report(
        args.checkv, args.abricate_plasmidfinder, args.trep_scaffold)
    PPdb.to_csv(f'{args.outbase}.sum.PhagePlasmid.csv', index=False)

    # 2) Create everything table
    stb = load_stb(args.stb)
    Edb = everything_report(
        stb, args.eukrep, PPdb, args.GTDB_table, args.checkm, args.trep_genome)
    Edb.to_csv(f'{args.outbase}.sum.AllGenomes.csv', index=False)
    assert len(set(Edb['domain']) - set(["Bacteriophage", "Plasmid", "Prokaryote", "Eukaryote"])) == 0

    # 3) Create everything stb
    for i, row in Edb[Edb["domain"].isin(["Bacteriophage", "Plasmid"])].iterrows():
        scaff = '_'.join(row['genome'].split("_")[1:])
        if scaff not in stb:
            stb[scaff] = row['genome'] + '.fa.gz'
    with open(f"{args.outbase}.sum.stb", 'w') as o:
        for s, b in stb.items():
            o.write(f"{s}\t{b}\n")

    # 3) Create structured dRep directory
    make_drep_directory(args.fasta, stb, Edb, args.outbase)

def make_drep_directory(fasta, stb, Edb, outbase):
    """
    Create structured dRep directory for  this sample
    """
    # 1) Create the base folder structure
    base_folder = f"dRep_sf_{outbase}"
    os.mkdir(base_folder)

    # Handle Euks without scoring
    Edb['contamination'] = [0 if ((c != c) & (d == 'Eukaryote')) else c for c, d in zip(Edb['contamination'], Edb['domain'])]
    Edb['completeness'] = [100 if ((c != c) & (d == 'Eukaryote')) else c for c, d in zip(Edb['completeness'], Edb['domain'])]

    # Handle plasmids without scoring
    Edb['contamination'] = [100 if ((c != c) & (d == 'Plasmid')) else c for c, d in zip(Edb['contamination'], Edb['domain'])]
    Edb['completeness'] = [0 if ((c != c) & (d == 'Plasmid')) else c for c, d in zip(Edb['completeness'], Edb['domain'])]

    for domains, name in zip([['Prokaryote', 'Eukaryote'], ['Bacteriophage', 'Plasmid']], ['set1', 'set2']):

        # 2) Create the folder
        s1f = os.path.join(base_folder, name + '/')
        os.mkdir(s1f)

        # 2.1) Create the genomeInfo.csv
        sdb = Edb[Edb['domain'].isin(domains)]
        sdb['gg'] = [g + '.fa.gz' for g in sdb['genome']]
        del sdb['genome']
        sdb = sdb.rename(columns={'gg':'genome', 'Contamination':'contamination', 'Completeness':'completeness'})

        sdb = sdb[['genome', 'contamination', 'completeness']]
        sdb = sdb.sort_values('completeness', ascending=False).drop_duplicates(subset=['genome'], keep='first')
        sdb.to_csv(os.path.join(s1f, f'{name}_genomeInfo.csv'), index=False)

        # 2.2) Create bins
        sbins = set(sdb['genome'])
        bf = os.path.join(s1f, 'bins/')
        os.mkdir(bf)
        stb_set1 = {s:b[:-6] for s, b in stb.items() if b in sbins}
        extract_bins(fasta, stb_set1, bf)
        cmd = f"gzip {bf}*.fa"
        subprocess.call(cmd, shell=True)

    # # 3) Create the set2 folder
    # s2f = os.path.join(base_folder, "set2")
    # os.mkdir(s2f)

    # # 3.1) Create the set2 genomeInfo.csv
    # sdb = Edb[Edb['domain'].isin(['Bacteriophage', 'Plasmid'])]
    # sdb['gg'] = [g + '.fa.gz' for g in sdb['genome']]
    # del sdb['genome']
    # sdb = sdb.rename(columns={'gg':'genome', 'Contamination':'contamination', 'Completeness':'completeness'})






def extract_bins(fasta, stb, out_base):
    if (out_base != '') & (out_base[-1] != '/'):
        out_base += '_'

    # Do it with append
    if fasta[-3:] == '.gz':
        handle = gzip.open(fasta, "rt")
        seqs = SeqIO.parse(handle, "fasta")

    else:
        seqs = SeqIO.parse(fasta, "fasta")

    for seq_record in seqs:
        id = str(seq_record.id).strip()
        if id not in stb:
            # print("{0} not in stb".format(id))
            continue
        fasta = stb[id]
        
        nw = open("{0}{1}.fa".format(out_base, fasta), 'a')
        nw.write('\n'.join([">{0}".format(id), str(seq_record.seq), '']))
        nw.close()

    if fasta[-3:] == '.gz':
        handle.close()

def phage_plasmid_report(checkv_folder, abr_loc, trep_loc):
    """
    All circular scaffolds are reported
    
    Phage are 'checkv_quality' in 'Medium-quality', 'High-quality', or 'Complete'
    
    Plasmids are cirulcar + at least one plasmid gene
    """
    # 1) Load a database of ciruclar contigs
    cir_loc = os.path.join(checkv_folder, 'complete_genomes.tsv') 
    Cdb = pd.read_csv(cir_loc, sep='\t')
    Cdb = Cdb[['contig_id', 'contig_length', 'prediction_type', 'repeat_length']].rename(columns={'contig_id':'scaffold', 'contig_length':'length', 'prediction_type':'repeat_type'})
    Cdb = Cdb[Cdb['repeat_type'] != 'Provirus'].sort_values('length', ascending=False).reset_index(drop=True)
    Cdb['circular'] = True
    
    # 2) Load checkv results
    checkv_loc = os.path.join(checkv_folder, 'quality_summary.tsv') 
    Vdb = pd.read_csv(checkv_loc, sep='\t')
    Vdb = Vdb[Vdb['checkv_quality'].isin(['Medium-quality', 'High-quality', 'Complete'])]
    Vdb = Vdb[['contig_id', 'contig_length', 'checkv_quality', 'miuvig_quality', 'completeness', 'contamination']].rename(columns={'contig_id':'scaffold', 'contig_length':'length'})
    Vdb['phage'] = True
    
    # 3) Load abricate results
    adb = pd.read_csv(abr_loc, sep='\t')
    adb = adb.groupby('SEQUENCE')['GENE'].agg('nunique').reset_index().sort_values('GENE').rename(columns={'SEQUENCE':'scaffold', 'GENE':'plasmid_genes'})
    adb['plasmid'] = True
    
    # 4) Incorporate taxonomy
    tdb = pd.read_csv(trep_loc, sep='\t')
    tdb = tdb[['scaffold', 'full_taxonomy', 'taxonomy']]
    
    # 5) Merge it all together
    CVdb = pd.merge(Cdb, Vdb, how='outer')
    CVAdb = pd.merge(CVdb, adb, how='left')
    Rdb = pd.merge(CVAdb, tdb, how='left')
    
    front = ['scaffold', 'length']
    metas = ['circular', 'phage', 'plasmid']
    for meta in metas:
        Rdb[meta].fillna(False, inplace=True)
    Rdb = Rdb[front + metas + [c for c in Rdb.columns if c not in set(front).union(set(metas))]].sort_values('length', ascending=False).reset_index(drop=True)
    
    return Rdb

def everything_report(stb, EukRep_table, phage_plasmid_table, GTDB_table, checkM_table, trep_table,
                     min_euk_p=50, min_pp_length=10000):
    """
    Create a full report on every binned genome
    
    De-replicate on the scaffold level (you might as well do this)
    
    Create a folder structure for the genomes, dereplicated on the scaffold level (this will be the input for the next dereplication step; could even create checkm / custom scoring inputs into the folder structure (with a folder like "input_tables"))
    """
    print("# 1) Load scaffold to bin from file")
    
    print("# 2) Load bacteria checkm info")
    bdb = pd.read_csv(checkM_table, sep='\t')
    bdb = bdb[(bdb['Completeness'] >= 50) & (bdb['Contamination'] < 10)]
    bdb = bdb[['Bin Id', 'Completeness', 'Contamination']].rename(columns={'Bin Id':'genome', 'Completeness':'completeness', 'Contamination':'contamination'})
    bdb['domain'] = 'Prokaryote'
    
    print("# 3) Load Euk info")
    try:
        edb = pd.read_csv(EukRep_table)
        edb = edb[edb['p_euk_length'] >= min_euk_p]
        edb['Bin Id'] = [str(c).split('.fa.gz')[0] for c in edb['bin']]
        edb['domain'] = 'Eukaryote'
        edb = edb[['Bin Id', 'domain', 'p_euk_length']].rename(columns={'Bin Id':'genome'})
    except:
        edb = pd.DataFrame()
    
    print("# 4) Add taxonomy for big genomes")
    Bdb = pd.concat([bdb, edb])
    
    print("# 4.1) Load gtbd")
    dbs = []
    for l in GTDB_table:
        dbs.append(pd.read_csv(l, sep='\t'))
    gdb = pd.concat(dbs).reset_index(drop=True)
    gdb = gdb[['user_genome', 'classification']].rename(columns={'user_genome':'genome', 'classification':'GTDB_classification'})
    gdb['GTDB_domain'] = [x.split(';')[0].replace('d__', '')for x in gdb['GTDB_classification']]
    Bdb = pd.merge(Bdb, gdb, how='left')
    
    print("# 4.2) Load tRep")
    tdb = pd.read_csv(trep_table, sep='\t')
    tdb['Bin Id'] = [str(c).split('.fa.gz')[0] for c in tdb['bin']]
    tdb = tdb[['Bin Id', 'full_taxonomy', 'taxonomy']].rename(columns={'Bin Id':'genome', 'taxonomy':'tRep_tax', 'full_taxonomy':'tRep_full_tax'})
    Bdb = pd.merge(Bdb, tdb, how='left')
    
    print("# 5) Load phage / plasmid info")
    if type(phage_plasmid_table) == type(pd.DataFrame()):
        pdb = phage_plasmid_table
    else:
        pdb = pd.read_csv(phage_plasmid_table)
    if len(pdb) > 0:
        pdb['domain'] = pdb.apply(call_domain, axis=1)
        pdb = pdb[~pdb['domain'].isna()]
        pdb['genome'] = [f"{d}_{s}" for d, s in zip(pdb['domain'], pdb['scaffold'])]
        pdb = pdb[pdb['length'] >= min_pp_length]
        pdb = pdb[['genome', 'full_taxonomy', 'taxonomy', 'domain', 'checkv_quality', 'circular', 'plasmid_genes', 'length', 'completeness', 'contamination']]
        pdb = pdb.rename(columns={'full_taxonomy':'tRep_full_tax', 'taxonomy':'tRep_tax', 'length':'scaffold_length'})
    
    print("# 6) Parse result table")
    Bdb = pd.concat([Bdb, pdb])
    front = ['genome', 'domain', 'tRep_full_tax']
    Bdb = Bdb[front + [c for c in Bdb.columns if c not in set(front)]].reset_index(drop=True)
    
    return Bdb

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

def load_stb_genomes(genomes):
    stb = {}
    for stb_file in genomes:
        bin_name = os.path.basename(stb_file).split('.')[1]
        if stb_file.endswith('.gz'):
            handle = gzip.open(stb_file, "rt")
        else:
            handle = open(stb_file, "r")

        for line in handle.readlines():
            if line[0] == '>':
                scaff = line.split()[0].strip()[1:]
                stb[scaff] = bin_name
    
    return stb

def compare_stb_files(stb, stb2):
    """
    Return a dicitonary of stb2 genome -> stb genome
    """
    db1 = pd.DataFrame(stb.items(), columns=['scaffold', 'bin1'])
    db2 = pd.DataFrame(stb2.items(), columns=['scaffold', 'bin2'])
    sdb = pd.merge(db1, db2)
    sdb = sdb[sdb['bin1'] != '0']
    for b1, db in sdb.groupby('bin1'):
        assert len(set(db['bin1'])) == 1
        if len(set(db['bin2'])) != 1:
            print('FAILURE!')
            print(db['bin2'].value_counts())
    return sdb

def call_domain(row):
    if row['phage']:
        return "Bacteriophage"
    elif (row['plasmid'] and row['circular']):
        return "Plasmid"
    else:
        return np.nan

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # General input
    parser.add_argument("-f", "--fasta",                    action="store",     help='Location of .fasta file binning was run on')
    parser.add_argument("-s", "--stb",                      action="store",     help='Location of stb file')
    parser.add_argument("-j", "--jgi_cov",                  action="store",     help="JGI coverage table")

    # Euk inputs
    parser.add_argument("-e", "--eukrep",                   action="store",     help='Location of EukRep output table')

    # Prokaryote inputs
    parser.add_argument("-g", "--GTDB_table",               action="store",     help='Location of GTDB tables', nargs='*')
    parser.add_argument("-chm", "--checkm",                 action="store",     help='Location of checkm table')

    # Phage / plasmid inputs
    parser.add_argument("-chv", "--checkv",                 action="store",     help='Location of checkv output folder')
    parser.add_argument("-a", "--abricate_plasmidfinder",   action="store",     help='Location of abricate plasmidfinder output')

    # Annotation inputs
    parser.add_argument("-ts", "--trep_scaffold",           action="store",     help='Location of tRep scaffold output')
    parser.add_argument("-tg", "--trep_genome",             action="store",     help='Location of tRep genome output')

    # Output
    parser.add_argument("-o", "--outbase",                      action="store",     help='Name of the output basename')

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
