
import pandas as pd
import os, sys
import glob



# Change to the directory with your shortstack output
os.chdir('/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack')

# Create path to (shortstack) Result files
path = '/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/*/Results.txt'

for result in glob.glob(path):
    # Creates the accession name, and the new file path in fname
    name = str(result)
    # Change the index in the list according to your path
    name = name[104:-12]
    fname = '/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/' + name + '/Results_miRNA.txt'
    
    # Creates df of result, add accession column, select MIRNA=Y and save in new file
    shs = pd.read_csv(result, sep="\t")
    Y = shs[shs['MIRNA'] == 'Y']
    Y.insert(0, 'Accession', name)
    Y.to_csv(fname, index=False, sep="\t")

# Create path to Result_miRNA files
path1 = '/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/*/*RNA.txt'

list_miRNA = [] 
for miRNA in glob.glob(path1):
    # Each miRNA file is transferred to a df, and subsequently appended to the list
    df = pd.read_csv(miRNA, sep="\t") 
    list_miRNA.append(df) 
# The lists are concatenated
concat = pd.concat(list_miRNA)

# The new df is saved in an output file
concat.to_csv('/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/concat_miRNA.txt', index=False, sep="\t")

