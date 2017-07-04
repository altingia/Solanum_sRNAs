This [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) pipeline is meant to plot results from the "srnadesc" Snakemake pipeline using [Circos](http://circos.ca/) .
*  It  will plot the sRNA-Seq coverage of sRNA on the reference genome chromosomes
*  It will plot the gene density and repeat density on the reference genome chromosomes 

# How to use it

## Install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html)

## Create a virtual environment named "myrnaseqworkflow" from the `environment.yaml` file
conda env create --name srnadesc2circosplot --file environment.yaml

## Activate this virtual environment
source activate srnadesc2circosplot

## Inputs needed (obtained preferentially from the srnadesc pipeline although not mandatory)
*  A tabulated text file containing the reference genome chromosome names and lengths
*  One GFF file for genes
*  One GFF file for repeats
*  One bedgraph file per sample  
*  One 

# Usage 
If Snakemake is properly set and if the pipeline is called `Snakefile` then just type `Snakemake` and it will run.
You can try a dry run (nothing happens, it checks the different steps and print the shell command) by typing `Snakemake -pn`

# Main outputs:
*  Circos plot (.png and .svg)Number of reads mapped at each step (before and after trimming, spikes, virus, rRNA, etc.)


 
