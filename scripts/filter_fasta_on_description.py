#!/usr/bin/env

"""
Filter repeat database to keep S. lycopersicum ribosomal RNA" Count the non-redundant sequence occurences in a fastq file

input: a fasta file
output: a filtered file

"""

####################
## Import librairies
####################
# import librairies
import argparse
from Bio import SeqIO

####################
## Arguments
####################
parser = argparse.ArgumentParser()
parser.add_argument("-i","--repeatfile",type=str,help="the fasta file containing repeat sequences to filter")
parser.add_argument("-o","--outfile",type=str,help="a filtered fasta file ")
parser.add_argument("-","--infilters",nargs="+",help="list of words to filter on",required=True)
args = parser.parse_args()


with open(args.repeatfile,"r") as filin:
	recs = [rec for rec in SeqIO.parse(filin,"fasta")]


# filter 
filtered_recs = []
for rec in recs:
	if any(word in rec.description for word in args.filters):
		filtered_recs.append(rec)

SeqIO.write(filtered_recs,args.outfile,"fasta")