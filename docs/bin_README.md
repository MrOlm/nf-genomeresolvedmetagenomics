# Overview of how bin works

1) Run prodigal
- This calls genes for all scaffolds in each assembly

2) Run bowtie2
- This does cross-mapping if available, or just maps reads to the assembly if not
- The only reason it does this is to provide coverage information to the binning algorithm

3) Run metabat2
- This is the primary binning algorithm being used. It groups scaffolds into bins

4) Run checkM on metabat2 bins
- This tries to estimate the completeness and contamination of all bins, assuming they're either archaea or bacteria

5) Run GTDB on metabat2 bins
- This tries to get the taxonomy of all bins, assuming they're either archaea or bacteria

6) Run bacteriophage binning
- This involves three programs:
    - virsorter2: tries to identify viral scaffolds
    - checkv: tries to estimate the completeness and contamination of all scaffolds
    - virfinder: tries to identify viral scaffolds

7) Run eukrep
- This tries to assess the eukaryotic precentage of each input scaffold AND each bin made by metabat2

8) Run abricate (to identify plasmids)
- This tries to find circular scaffolds, which can be plasmids

9) Preform gene annotation
- This involves two programs:
    - DIAMOND: aligns amino acid sequences to the UniRef databases
    - tRep: takes those DIAMOND alignments and parses the results

10) Summarize binning
- This is the custom program "summarize_binning.py" that is meant to take in all of the above information and summarize it in a helpful way

# Details of specific aspects

## Sumarize binning

This python program does a couple of things and feel free to check it out yourself (it's in the `bin` folder of this repository)

1) Create phage / plasmid table "sum.PhagePlasmid.csv"
    - 1) Load a database of ciruclar contigs from checkV (complete_genomes.tsv)
    - 2) Load checkV results (quality_summary.tsv). Things assigned as 'Medium-quality', 'High-quality', and 'Complete' are considered phage
    - 3) Load abricate results. All scaffolds with at least one plasmid gene are considered plasmids
    - 4) Incorporate taxonomy. This comes straight from  trep
    - 5) Merge it all together into one table

2) Create everything table ".sum.AllGenomes.csv"
    - 1) Bacteria is anything with "['Completeness'] >= 50) & (bdb['Contamination'] < 10" according to checkM
    - 2) Euks are anything with p_euk_length > 50
    - 3) Add GTDB taxonomy (called GTDB_classification)
    - 4) Add tRep taxonomy (called tRep_tax)
    - 5) Merge in the PhagePlasmid table referenced above

3) Do some scaffold-level-dereplication
    - If a scaffold is called a phage but is also in a bacterial genome, no longer call it a phage (for example)
    - Right now this impacts "sum.stb"

4) Make a structured directory for dRep
    - Make some folders to run dRep on based on the above information