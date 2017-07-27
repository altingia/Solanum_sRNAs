"""
read_length_distribution: returns the length of sequences from a fasta or fastq file



"""

# import libraries
from Bio import SeqIO
import argparse

####################
## Arguments
####################
parser = argparse.ArgumentParser()
parser.add_argument("-f","--inputfile",type=str,help="the fasta or fastq file from where sequence lenght will be assessed")
parser.add_argument("-o","--outfile",type=str,help="the tabulated text file returned ",default="nr_seq_counts.txt")
parser.add_argument("-n","--name",type=str,help="name of the sample (will be used as the column name",default="sample")
args = parser.parse_args()


if args.inputfile.endswith((".fasta",".fa")):
	with open(args.inputfile,"r") as filin:
		seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fasta")]
if args.inputfile.endswith((".fastq",".fq")):
	with open(args.inputfile,"r") as filin:
		seqs = [str(rec.seq) for rec in SeqIO.parse(filin,"fastq")]

# measure sequence length
with open(args.outfile,"w") as fileout:
	fileout.write(args.name + "\n")
	for seq in seqs:
		fileout.write(str(len(seq)) + "\n")



