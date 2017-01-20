# !/usr/bin/env python
"""
To count the number of reads at a certain stage in the pipeline, I created a Python script that counts the number of sequences.
It returns a simple line separated by tabulation:
"sample name"	\t	"Step name e.g. 'after trimming'"	\t "number of counts"
"LA1777"	\t	"After trimming"			\t "127826"			

This Python function reads a fastq file and returns a text file containing a single line containing:
* --field1		a name that can be specified by the user (default = fastq file name)
* --field2 Field n°2: the stage of interest 
* Field n°3: the number of fastq records (equal to sequences)
"""

# import modules
import argparse
from Bio import SeqIO

# define a parser 
parser = argparse.ArgumentParser(description="parser for counting fastq records")

parser.add_argument("-f","--fastqfile",type=str,help="the fastq file to count from")
parser.add_argument("-o","--outfile",type=str,help="the tabulated text file to write results")
parser.add_argument("-f1","--field1",type=str,help="what to write in the first field e.g. file name",default="sample")
parser.add_argument("-f2","--field2",type=str,help="what to write in the second field e.g. analysis step",default="step")
args = parser.parse_args()

# function
def count_fastq_seqs(infile,outfile,sample,stepname):
    with open(args.fastqfile,"r") as filin:
        seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fastq")]
        nseqs = len(seqs)
    with open(args.outfile,"w") as fileout:
        fileout.write(str(sample) + "\t" + str(stepname) + "\t" + str(nseqs) + "\n")   

count_fastq_seqs(args.fastqfile,args.outfile,args.field1,args.field2)
