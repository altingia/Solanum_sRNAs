import pandas as pd
import os
import sys
import glob

class SHS(object):
    def __init__(self, shs_i_path, path, column, value):
        self.shs_i_path = shs_i_path
        self.path = path
        self.column = column
        self.value = value
        
    def select_data(self):
        # from source file a column and certain value can be selected. A column with the source accession is added. For each accession a file is created.
        path = self.shs_i_path
        os.makedirs(self.path + '/' + self.column + '_' + self.value)
        for result in glob.glob(path):
            df = pd.read_csv(result, sep="\t")
            selection = df[df[self.column] == self.value]
            head = os.path.dirname(result)
            name = os.path.basename(head)
            selection.insert(0, 'Accession', name)
            fname = self.path + '/' + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + name + '.txt'
            selection.to_csv(fname, index=False, sep="\t")
    
    def concat_data(self):
        # to concat files
        path = self.path + '/' + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + '*' + '.txt'
        list_i=[]
        for i in glob.glob(path):
            df = pd.read_csv(i, sep="\t") 
            list_i.append(df)
        concat = pd.concat(list_i)
   
        o_path = os.path.dirname(path) + "/" + self.column + "_" + self.value + "_concat.txt"
        concat.to_csv(o_path, index=False, sep="\t")
    
    def pres_in_accs(self):
        # to add a column in which the presence of a certain miRNA in accessions is decribed
        path = self.path + '/' + self.column + '_' + self.value + '/' + self.column + '_' + self.value + '_' + 'concat.txt'
        df = pd.read_csv(path, sep="\t")
        df.insert(1, 'Present in accessions', 'unique')
        uni = df.drop_duplicates(['MajorRNA'], keep=False)
        red = df[df.duplicated(['MajorRNA'], keep=False)]

        RNA_list = []
        for mol in red['MajorRNA']:
            if mol not in RNA_list:
                RNA_list.append(mol)
            else:
                pass
        list_red = []
        for mol in RNA_list:
            list_red.append(red[red['MajorRNA'] == mol])

        length_list_red = len(list_red)

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
            dict_RNA_in_accs.update(lib)

        for i in range(len(df)):
            current_RNA = str(df['MajorRNA'][i])
            if current_RNA in dict_RNA_in_accs:
                df.set_value(i, 'Present in accessions', dict_RNA_in_accs[current_RNA])
            else:
                df.set_value(i, 'Present in accessions', df['Accession'][i])

        o_path = os.path.dirname(path) + "/" + self.column + "_" + self.value + "_pres_in_accs.txt"
        df.to_csv(o_path, index=False, sep="\t")
