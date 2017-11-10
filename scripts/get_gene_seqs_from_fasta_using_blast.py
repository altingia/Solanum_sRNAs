"""
Goal = using a single fasta "bait" sequence and a RNA-Seq de novo assembly (multifasta) -> collects the best Blast hit and writes it to a fasta file
1) make a blastdb for the de novo RNA-Seq assembly
2) performs a blast using the bait and assembly
2) filter the output (keep best hit ID)
3) retrieve the corresponding fasta file
4) write to a fasta file

[Usage]
python get_seqs -b baits/ -a assemblies/ -outdir results/ -blast blast_results/
"""

import os
import subprocess
import argparse
import fnmatch


# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b","--baits",type=str,help="a directory containing fasta files to search in the assemblies",default="baits/")
parser.add_argument("-a","--assemblies",type=str,help="a directory containing RNA-Seq de novo assemblies",default="assemblies/")
parser.add_argument("-o","--outdir",type=str,help="a directory to place the extracted fasta results",default="results/")
parser.add_argument("-blast","--blastresdir",type=str,help="a directory to place the blast results",default="blast_results/")
args = parser.parse_args()

# Collect bait genes and assemblies
assemblies = [assembly for assembly in os.listdir(args.assemblies) if assembly != ".DS_Store" and assembly.endswith(".fasta")]
genes = [gene for gene in os.listdir(args.baits) if gene != ".DS_Store"]


###########################################################
# Creates result and blast directories if they don't exist
##########################################################
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

if not os.path.exists(args.blastresdir):
    os.makedirs(args.blastresdir)

####################################
# for each assembly, makes a blastdb
####################################
for assembly in assemblies:
	if os.path.exists(args.assemblies + assembly + ".nhr"):
		pass
	else:
		# make database
		blastdb_command = "makeblastdb -in " + args.assemblies + assembly + " -dbtype nucl"
		subprocess.call(blastdb_command,shell=True)

######################
# index the fasta file
######################
for assembly in assemblies:
	samtools_index = "samtools faidx " + args.assemblies + assembly
	subprocess.call(samtools_index,shell=True)

###############
# Get bait fasta, de novo RNA-Seq assembly and blastn
###############
#get bait and assembly
for bait in genes:
	for assembly in assemblies:
		accession = os.path.splitext(assembly)[0]
		gene = os.path.splitext(bait)[0]
		# perform blast
		blastresfile = args.blastresdir + accession + "_" + gene + ".outfmt6"
		blast_command = "blastn -db " + args.assemblies + assembly + " -query " + args.baits + bait + " -outfmt 6 -out " + blastresfile 
		subprocess.call(blast_command,shell=True)

		outfilename = args.outdir + accession + "_" + gene + ".fasta" 
		#check if the file is empty		
		#if os.stat(blastresfile).st_size == 0:		
#			with open(outfilename,"w") as fileout:
#				fileout.write("no hit found")
#		else:
		#filter the output (keep best hit ID)
		with open(blastresfile,'r') as infile:
			first_line = infile.readline()
			items = first_line.split("\t")
			seq2retrieve = items[1] # sequence id
			# retrieve corresponding sequence in the RNA-Seq de novo assembly and write to fasta
			samtools_retrieve = "samtools faidx " + args.assemblies + assembly + " " + seq2retrieve + " > " + outfilename
			subprocess.call(samtools_retrieve,shell=True)

#########################################
## Creates a compiled multifasta per gene
#########################################
for gene in genes:
	fastas = [args.outdir + f for f in os.listdir(args.outdir) if gene in f]
	multifasta_for_one_gene = args.outdir + gene
	with open(multifasta_for_one_gene,"w") as outfile:
		for fasta in fastas:
			with open(fasta) as infile:
				for line in infile:
					outfile.write(line)