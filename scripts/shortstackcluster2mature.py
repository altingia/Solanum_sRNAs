# /usr/bin/env python 


"""
Takes a ShortStack cluster and extracts the mature RNA sequence + sends it to a fasta file
[Usage] python shortstackcluster2mature.py -i [cluster text file] -o [mature fasta file]
"""

##################
# import libraries
##################
import argparse

###########
# arguments
###########

parser = argparse.ArgumentParser()
parser.add_argument("-i","--clusterfile",type=str,help="the ShortStack MIRNA cluster text file to convert")
parser.add_argument("-o","--mature_fasta",type=str,help="the output fasta file for the mature miRNA from the selected cluster")
args = parser.parse_args()

with open(args.clusterfile,"r") as filin, open(args.mature_fasta,"w") as fileout:
	for i, line in enumerate(filin):
		if i==0:
			# write the name of the cluster
			fileout.write(">" + line.strip().split(" ")[0] + "\n")
		if i==4:
			if "miRNA " in line:
				# removes the dots
				mirna_seq = line.strip().replace(".","")	
				# split string into a list to get the first element (always the miRNA)
				l = mirna_seq.split(" ")
				mirna = l[0]
				# print to a fasta file
				fileout.write(mirna + "\n")