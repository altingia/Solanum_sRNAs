#!/usr/bin/env

"""
Count the non-redundant sequence occurences in a fastq file

input: a fastq file
output: a tabulated file containing two columns

Example:
>read1
ATCGATCG
>read2
ATCGATCG
>read3
ATCGATCG

Then "ATCGATCG" will be counted 3 times

Column1	Column2
ATCGATCG	3
ATCGAAAA	10
TCGTGCTG	1000
"""

####################
## Import librairies
####################
# import librairies
import pandas as pd
import argparse
import collections
from collections import Counter
from Bio import SeqIO

####################
## Arguments
####################
parser = argparse.ArgumentParser()
parser.add_argument("-i","--fastq",type=str,help="the fastq file from where sequences will be counted")
parser.add_argument("-o","--outfile",type=str,help="the tabulated text file returned ",default="nr_seq_counts.txt")
args = parser.parse_args()


with open(args.fastq,"r") as filin:
            seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fastq")]
            c = Counter(seqs)
            df = pd.DataFrame(list(c.items()),columns=["seq","counts"])
            #df.sort(["counts"],ascending=[False],inplace=True)
            df.sort_values(by="counts",ascending=False,inplace=True)
df.to_csv(args.outfile,sep="\t",index=False) 