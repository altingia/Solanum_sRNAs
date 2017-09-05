# /usr/bin/env python 


"""
For one sample, from the genes predicted by sPARTA to be targets, gets their corresponding BED file of their coordinates

Usage:
python sparta_targets_to_bed_file.py -i All.libs.validated.uniq.csv -g ITAG3.0_gene_models.gff -o genes.bed
"""

##################
# import libraries
##################
import argparse
import pandas as pd

###########
# arguments
###########

parser = argparse.ArgumentParser()
parser.add_argument("-i","--sparta",type=str,help="the sPARTA All.libs.validated.uniq.csv file to be converted")
parser.add_argument("-g","--gff",type=str,help="the GFF annotation file from where the coordinates will be taken")
parser.add_argument("-o","--out",type=str,help="the BED file parsed",default="genes.bed")
args = parser.parse_args()

# Get the genes from sPARTA 
sparta = pd.read_csv(args.sparta)
target_genes= sparta["Target"].tolist()

# Import the GFF file
# Keep gene features only
# Creates a new column with the gene ID
df = pd.read_csv(args.gff,comment='#',sep="\t",header=None)
df.columns = ["chr","source","feature","start","end","score","strand","frame","attribute"]
df = df[df['attribute'].str.contains("ID=gene:")]

#df = df[df.feature == "gene"]
attributes = df["attribute"].tolist()
attributes = [att.split(";")[0].replace("ID=","") for att in attributes]
df["gene_id"] = attributes

# Filters the GFF file to keep target_genes
# Write to disk
df = df[df["gene_id"].isin(target_genes)]

# Convert to bed
# 1. chrom
# 2. chromStart
# 3. chromEnd
# 4. name 
# 5. score
# 6. strand
# to bed
bed = df[["chr","start","end","gene_id","score","strand"]]
new_starts = bed["start"] - 1
bed1 = bed.copy()
bed1["new_starts"] = new_starts
bed1 = bed1.drop("start",1) # delete old start positions
bed1.head()
new_col_ordering = ["chr","new_starts","end","gene_id","score","strand"]
bed = bed1[new_col_ordering]
bed.head()
bed.to_csv("parsed.gff",sep="\t",index=False,header=False)
