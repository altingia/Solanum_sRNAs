#!/usr/bin/env python

"""
Take one Shortstack result file and creates a fasta file containing all miRNAs
Usage:
python get_mirnas_from_shortstack_res.py -i [shortstack result file] -o [path/to/outfile] 

Example: python get_mirnas_from_shortstack_res.py -i shortstack/C32/Results.txt -o C32_miRNAs.fasta 
"""

####################
## Import librairies
####################
# import librairies
import pandas as pd
import argparse

####################
## Arguments
####################
parser = argparse.ArgumentParser()
parser.add_argument("-i","--shortstackfile",type=str,help="the Shortstack result file")
parser.add_argument("-o","--outfastafile",type=str,help="the fasta file containing miRNAs sequences returned as an output",default="miRNAs.fasta")
args = parser.parse_args()

#########################################
## Import all ShortStack results
##########################################

# read and store each ShortStack result file as a pandas dataframe
df = pd.read_table(args.shortstackfile)

# select only miRNAs
df = df[df.MIRNA == "Y"]
df = df[["Name","MajorRNA"]]

###########################
## Get non-redundant miRNAs
###########################
# convert to list
miRNAs = list(df["MajorRNA"].tolist())
names = list(df["Name"].tolist())

# write to file
with open(args.outfastafile,"w") as fileout:
	for i in range(0,len(miRNAs),1):
		fileout.write(">" + names[i] + "\n" + miRNAs[i] + "\n")

