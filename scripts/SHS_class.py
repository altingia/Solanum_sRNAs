import pandas as pd
import os
import sys
import glob

class SHS(object):
    def __init__(self, shs_i_path, path, column, value): 
                 #*args, **kwargs):
        #column = kwargs.get('column', None)
        #value = kwargs.get('value', None)
        self.shs_i_path = shs_i_path
        self.path = path
        self.column = column
        self.value = value
        
def select_data(self):
    # from source file a column and certain value can be selected. A column with the source accession is added. 
    # For each accession a file is created.
    path = self.shs_i_path
    
    if self.column == '':
        if os.path.exists(self.path + '/all') == False:
            os.makedirs(self.path + '/all')
        else:
            print("directory already exists")
            
        for result in glob.glob(path):
            df = pd.read_csv(result, sep="\t")
            head = os.path.dirname(result)
            name = os.path.basename(head)
            fname = self.path + '/all/Results_' + name + '.txt'
            df.to_csv(fname, index=False, sep="\t")
    else:
        if os.path.exists(self.path + '/' + self.column + '_' + self.value) == True:
            print("directory already exists")
        else:
            os.makedirs(self.path + '/' + self.column + '_' + self.value)
            print("directory is created")
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
        print(mol)
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
        print(RNA_red)
        RNA = ""
        for rna in RNA_red:
            if rna not in RNA:
                RNA = rna
        acc_red = list_red[i]['Accession']
        print(acc_red)
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


def unique(self, accession):
    path = self.path + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + accession + '.txt'
    df = pd.read_csv(path, sep="\t")
    uni = df.drop_duplicates(['MajorRNA'], keep=False)
    o_path = os.path.dirname(path) + "/" + self.column + "_" + self.value + "_uni_" + accession + ".txt"
    uni.to_csv(o_path, index=False, sep="\t")


def redundant(self, accession):
    path = self.path + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + accession + '.txt'
    df = pd.read_csv(path, sep="\t")
    red = df[df.duplicated(['MajorRNA'], keep='first')]
    o_path = os.path.dirname(path) + "/" + self.column + "_" + self.value + "_red_" + accession + ".txt"
    red.to_csv(o_path, index=False, sep="\t")

def bedfile(self):
    #path = self.path + '/' + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + '*' + '.txt'
    #dirname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(self.shs_i_path)))) 
    #diname = os.path.basename(dirname)
    #dname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(dirname))))
    #print(path)
    #print(dirname)
    #print(diname)
    #print(dname)
    if self.column == '':
        path = self.path + '/all/Results_' + '*' + '.txt'
        dirname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(self.shs_i_path)))) 
        diname = os.path.basename(dirname)
        dname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(dirname))))
        list_i=[]
        if os.path.exists(dname + '/11.analysis/sRNA/' + diname) == False:
            os.makedirs(dname + '/11.analysis/sRNA/' + diname)
            print("directory is created")
        else:
            print("head directory already present")
       
        if os.path.exists(dname + '/11.analysis/sRNA/' + diname + '/all')==False:
             os.makedirs(dname + '/11.analysis/sRNA/' + diname + '/all')
        else:
            print("sub directory already present")

        for i in glob.glob(path):
            bed = pd.read_csv(i, sep ='\t')           
            s1=bed['#Locus'].apply(lambda x: x.split(':'))
            bed['Chr'] = s1.apply(lambda x: x[0])
            bed['Range'] = s1.apply(lambda x: x[1])
            s2=bed['Range'].apply(lambda x: x.split('-'))
            bed['Start'] = s2.apply(lambda x: x[0])
            bed['End'] = s2.apply(lambda x: x[1])
            bed= bed.drop('#Locus', 1)
            bed= bed.drop('Range', 1)
            bedF= bed[['Chr', 'Start', 'End', 'Name', 'Strand']]
        
            name = os.path.basename(i)
            name = name.replace(".txt","")
            fname = dname + '/11.analysis/sRNA/' + diname + '/all/' + name + '.bed'
            print(fname)
            bedF.to_csv(fname, index=False, sep="\t")
            list_i.append(bed)
    
    else:
        path = self.path + '/' + self.column + '_' + self.value + '/Results_' + self.column + '_' + self.value + '_' + '*' + '.txt'
        dirname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(self.shs_i_path)))) 
        diname = os.path.basename(dirname)
        dname = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(dirname))))
        list_i=[]
        if os.path.exists(dname + '/11.analysis/sRNA/' + diname) == False:
            os.makedirs(dname + '/11.analysis/sRNA/' + diname)
            print("directory is created")
        else:
            print("head directory already present")
       
        if os.path.exists(dname + '/11.analysis/sRNA/' + diname + '/' + self.column + '_' + self.value )==False:
             os.makedirs(dname + '/11.analysis/sRNA/' + diname + '/' + self.column + '_' + self.value )
        else:
            print("sub directory already present")

        for i in glob.glob(path):
            bed = pd.read_csv(i, sep ='\t')           
            s1=bed['#Locus'].apply(lambda x: x.split(':'))
            bed['Chr'] = s1.apply(lambda x: x[0])
            bed['Range'] = s1.apply(lambda x: x[1])
            s2=bed['Range'].apply(lambda x: x.split('-'))
            bed['Start'] = s2.apply(lambda x: x[0])
            bed['End'] = s2.apply(lambda x: x[1])
            bed= bed.drop('#Locus', 1)
            bed= bed.drop('Range', 1)
            bedF= bed[['Chr', 'Start', 'End', 'Name', 'Strand']]
        
            name = os.path.basename(i)
            name = name.replace(".txt","")
            fname = dname + '/11.analysis/sRNA/' + diname + '/' + self.column + '_' + self.value + '/' + name + '.bed'
            print(fname)
            bedF.to_csv(fname, index=False, sep="\t")
            list_i.append(bed)
