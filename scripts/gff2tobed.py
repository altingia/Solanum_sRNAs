#!/usr/bin/env python

"""
Convert GFF2 file format to BED
Careful: GFF2 files have 1-based coordinates and BED files have 0-based coordinates
This is called "off-by-one coordinate"
"""

####################
## Import librairies
####################
# import librairies
import argparse
import os

####################
## Arguments
####################
parser = argparse.ArgumentParser()
parser.add_argument("-i","--gff2",type=str,help="the input GFF2 file to be converted")
parser.add_argument("-o","--bed",type=str,help="the output BED file returned",default="bed.txt")
args = parser.parse_args()



############

with open(args.gff2,"r") as filin:
    lines = filin.readlines()
    lines = [line.strip() for line in lines]
    newlines = [line.split("\t") for line in lines]
    
with open(args.bed,"w") as fileout:
    for i in range(0,len(lines),1):
        chrom = newlines[i][0]
        start = newlines[i][3]
        end = newlines[i][4] 
        name = newlines[i][8]
        score = newlines[i][5]
        strand = newlines[i][6]
        fileout.write(str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(name) + "\t" + str(score) + "\t" + str(strand) + "\n")


