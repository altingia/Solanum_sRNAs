#!/usr/bin/env python

"""
Merge all miRNA from Shortstack: we get a concatenated dataframe and a fasta file containing all miRNAs
Usage:
python merge_mirnas_from_shortstack.py -d [list of shortstack result files] -o [path/to/outfile] 

file hierarchy
shortstack/
	C32/
		|-- Results.txt
	LA4024/
		|-- Results.txt
	PI127826/
		|-- Results.txt

Example: python merge_mirnas_from_shortstack.py -d shortstack/C32/Results.txt shortstack/LA4024/Results.txt 
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
parser.add_argument("-l","--shortstackfiles",nargs="+",help="the Shortstack result files stored as a list")
parser.add_argument("-o","--outtextfile",type=str,help="the tabulated text file returned as an output",default="miRNAs.txt")
parser.add_argument("-f","--outfastafile",type=str,help="the fasta file containing miRNAs sequences returned as an output",default="miRNAs.fasta")
args = parser.parse_args()

#########################################
## Import all ShortStack results
##########################################

# read and store each ShortStack result file as a pandas dataframe
dfs = [pd.read_table(f) for f in args.shortstackfiles]

# select only miRNAs
for i in range(0,len(dfs),1):
	dfs[i] = dfs[i][dfs[i].MIRNA == "Y"]
	dfs[i] = dfs[i][["MajorRNA"]]

print("#############################")
print("miRNA column was selected. This is how one dataframe looks like")
print(dfs[0].head())
print("We have {0} miRNA files".format(len(dfs)))
print("##############################")

#########################################
## Start recursive merging (full outer join)
##########################################

# initialize the first dataframe and merge
df = pd.merge(dfs[0],dfs[1],on="MajorRNA",how="outer")

# starts to iterate and recursively merge the dfs
for i in range(2,len(dfs)-1,1):
    df = pd.merge(df,dfs[i],on="MajorRNA",how="outer")
    print("\n")
    print("#############################")
    rows, columns = df.shape
    print("This is the {0}th merge. The corresponding dataframe has {1} lines".format(str(i),str(rows)))
    print("##############################")

# the final iteration is to merge the last dataframe
df = pd.merge(df,dfs[len(dfs)-1],on="MajorRNA",how="outer")   
print("\n")
print("#############################")
rows, columns = df.shape
print("This is the {0}th merge. The corresponding dataframe has {1} lines".format(str(len(dfs)-1),str(rows)))
print("##############################")
print("\n")
print("merging done")


###########################
## Get non-redundant miRNAs
###########################
# convert to list
miRNAs = list(set(df["MajorRNA"].tolist()))
print("We found {0} non-redundant miRNAs".format(len(miRNAs)))

# write to file
with open("miRNAs.txt","w") as fileout1:
	fileout1.write("id"+"\t"+"miRNA_seq" + "\n")
	for i in range(0,len(miRNAs),1):
		fileout1.write(str(i+1)+"\t"+miRNAs[i]+"\n")
with open("miRNAs.fasta","w") as fileout2:
	for i in range(0,len(miRNAs),1):
		fileout2.write(">miRNA"+str(i+1)+"\n"+miRNAs[i]+"\n")

