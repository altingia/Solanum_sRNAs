import pandas as pd
import os, sys
import glob

# dataframe of input file
df = pd.read_csv('/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/concat_miRNA.txt', sep="\t")

# An extra column named "Present in accessions" with the value "unique" is created next to Accessions column
df.insert(1, 'Present in accessions', 'unique')

# Creates a df with ONLY the unique MajorRNA
uni = df.drop_duplicates(['MajorRNA'], keep=False)
# Creates a df with ONLY redundant MajorRNA
red = df[df.duplicated(['MajorRNA'], keep=False)]

# Creates a list of the different MajorRNA present in the red df. So each RNA sequence is listed only ones.
RNA_list = []
for mol in red['MajorRNA']:
    if mol not in RNA_list:
        RNA_list.append(mol)
    else:
        pass

# Then iterates over this RNA_list to make a list_red in which the accessions in which the concerning MajorRNA is present, are given.
list_red = []
for mol in RNA_list:
    list_red.append(red[red['MajorRNA'] == mol])

length_list_red = len(list_red)

# Then a dict is created with MajorRNA as key and accessions in which present as values
dict_RNA_in_accs = {}
for i in range(length_list_red):
    RNA_red = list_red[i]['MajorRNA']
    RNA = ""
    for rna in RNA_red:
        if rna not in RNA:
            RNA = rna
    acc_red = list_red[i]['Accession']
    list_acc_red = []
    for acc in acc_red:
        if acc not in list_acc_red:
            list_acc_red.append(acc)
    lib = {RNA: list_acc_red}
    dict_RNA_in_accs.update(lib)   #dict_RNA_in_accs is a dict of dicts

# This replaces the default value of column "Present in accessions" with the list of all accession in which the MajorRNA is present.
for i in range(len(df)):
    current_RNA = str(df['MajorRNA'][i])
    if current_RNA in dict_RNA_in_accs:
        df.set_value(i, 'Present in accessions', dict_RNA_in_accs[current_RNA])
    else:
        df.set_value(i, 'Present in accessions', df['Accession'][i])

df.to_csv('/Users/michelle/surfdrive/Shared/trichome_team/06.results_from_xp/sRNA-Seq/20170117_srnadesc/shortstack/summary_miRNA.txt', index=False, sep="\t")




p
