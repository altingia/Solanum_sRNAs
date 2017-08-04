# /usr/bin/env python 


"""
Takes a ShortStack cluster and extracts the hairpin RNA + sends it to a fasta file
[Usage] python shortstackcluster2hairpin.py -i [cluster text file] -o [hairpin fasta file]
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
parser.add_argument("-o","--hairpin_fasta",type=str,help="the output fasta file for the hairpin precursor of the miRNA from the selected cluster")
args = parser.parse_args()

with open(args.clusterfile,"r") as filin, open(args.hairpin_fasta,"w") as fileout:
	for i, line in enumerate(filin):
		if i==0:
			# write the name of the cluster
			fileout.write(">" + line.strip().split(" ")[0] + "\n")
		if i==2:
			fileout.write(line.strip() + "\n")

