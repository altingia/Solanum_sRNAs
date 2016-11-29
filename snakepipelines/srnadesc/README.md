This [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) pipeline is meant to characterize in details small RNA sequencing files.
*  It  will generate statistics on the number of small RNAs that aligns to different databases (phytoviruses, RFAM, repeats, transposons, etc.)
*  It will characterize the type of small RNAs that can be found in the samples (microRNAs, siRNAs, phased siRNAs, etc.) using specialized programs like [ShortStack](https://github.com/MikeAxtell/ShortStack) 

# Installation
## Snakemake
Please visit the [Snakemake homepage](https://bitbucket.org/snakemake/snakemake/wiki/Home) for detailed instructions about Snakemake.

## Dependencies
*  Bowtie (v1 for ShortStack) and Bowtie2 (for other aligments)
*  samtools v1.2.1
*  RNAfold from Vienna 

# Usage 
If Snakemake is properly set and if the pipeline is called `Snakefile` then just type `Snakemake` and it will run.
You can try a dry run (nothing happens, it just checks the different steps) by typing `Snakemake -n`

 
