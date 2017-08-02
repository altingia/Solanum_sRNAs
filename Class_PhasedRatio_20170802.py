import pysam
import sys
import os
import pandas as pd
import numpy as np
import glob
import math

class PhasedRatio(object):
    def __init__(self, sam_path, bed_path, shs_path, m): 
        self.sam_path = sam_path
        self.bed_path = bed_path
        self.shs_path = shs_path
        self.m = m

def abund_sRNA(self):
    path = self.sam_path
    ab= {}
    for result in glob.iglob(path):
        df = pd.read_csv(result, sep="\t", header = None, 
                 names = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
                          'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'TAG', 'TAG1', 'TAG2', 'TAG3', 'TAG4', 'TAG5', 'TAG6'])
    
        fname = os.path.basename(result)
        fname = fname.replace('nhC32_', '')
        fname = fname.replace('.sam', '')
        
        ab_all = len(df.FLAG)
        
        list_seq = []
        for seq in df.SEQ:
            if seq not in list_seq:
                list_seq.append(seq)
        ab_uni = len(list_seq)
        
        ab[fname] = int(ab_all), int(ab_uni)
    return(ab)

def phased(self):
    path = self.sam_path
    path2= self.bed_path
    m = int(self.m)
    
    directory = os.path.dirname(os.path.dirname(os.path.dirname(path))) + '/Bin/Cluster_Bin' + str(m)
    if os.path.exists(directory) == False:
        os.makedirs(directory)
    else:
        print("directory already exists")
    

    for result in glob.iglob(path):
        # For each samfile the cluster name is taken. In the bedfile the index corresponding
        # to the cluster is searched and then the start position of the cluster is determined.
        fname = os.path.basename(result)
        fname = fname.replace('nhC32_', '')
        fname = fname.replace('.sam', '')
        print(fname)
        
        start = pd.read_csv(path2, sep="\t")
        starti = np.where(start.Name == fname)
        print(starti)
        starti = str(starti)
        print(starti)
        starti = starti.replace('(array([', '')
        print(starti)
        starti = starti.replace(']),)', '')
        print('starti' + str(starti))
        spos = start.columns.get_loc('Start')
        print('spod' + str(spos))
        start = start.ix[int(starti),spos]

        # With the samfile for each read the '5 position is used to calculate the bin.
        # then the column with the bin/read is added to the df and saved in a file.
        DF = pd.read_csv(result, sep = '\t', header = None, 
                     names = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
                              'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'TAG', 
                              'TAG1', 'TAG2', 'TAG3', 'TAG4', 'TAG5', 'TAG6'])
    
    
        # optional to filter on CIGAR, the in phasing with m only reads with length m are included
        df = DF.copy()
        #cigar_filter = str(m) + 'M'
        #df = df[df.CIGAR == cigar_filter]
        #df = df[df.CIGAR == '24M']
        #df = df.reset_index(drop=True)
        #num = len(df.FLAG)
        #if num <= 4:
        #    continue
        
        
        numb = len(df.FLAG)
        insert_pos = df.columns.get_loc("POS")
        bin_list=[]
        for n in range(numb):
            numb = numb -1
            if df.ix[int(numb),1] == 0:
                pos = df.ix[int(numb),insert_pos]
                i = pos-(start-1)
                Bin = i % m
                bin_list.append(Bin)
            else:
                pos = df.ix[int(numb),insert_pos]
                i = pos-(start-1)
                if i-18 <= 0:
                    ii = 21 + (i-18)
                else:
                    ii = i-18
                Bin = ii % m
                bin_list.append(Bin)
        for n,i in enumerate(bin_list):
            if i==0:
                bin_list[n]=21
                
        rev_bin_list = bin_list[::-1]
        df.insert(insert_pos, 'BIN', rev_bin_list)
        out1name = str(fname + '.txt')
        df.to_csv(out1name, index=False, sep="\t")
        
        # For each bin the amount of reads within the bin is calculated.
        # Subsequently the total amount of reads is calculated and the phased ratio.
        bindf = pd.DataFrame({'BIN': range(1,22)})
        numb = len(bindf)+1
        bin_n = []
        for n in reversed(range(numb)):
            if n >= 1:
                amount = sum(df.BIN == n)
                bin_n.append(amount)
        
        rev_bin_n = bin_n[::-1]
        bindf['Reads'] = rev_bin_n
        # add column with DistinctReads
        dist = []
        for n in reversed(range(numb)):
            if n >= 1:
                temp_df = df[df.BIN == n]
                list_seq = []
                for seq in temp_df.SEQ:
                    if seq not in list_seq:
                        list_seq.append(seq)
                ab_uni = len(list_seq)
                dist.append(ab_uni)
        rev_dist = dist[::-1]
        bindf['DistinctReads'] = rev_dist
        

        #with this script the PR per highest or two highest bins is calculated. Column with highest bin Y/N is added. 
        #DF are pivotd and the total and PR are added.
        total = bindf['Reads'].sum()
        high = list(bindf['Reads'].nlargest(n = 3, keep = 'last'))
        multiple_max= []
        for n in high:
            if n in multiple_max:
                multiple_max.append('d')
            elif n not in multiple_max:
                multiple_max.append(n)

        length = len(multiple_max)
        if length == 3:
            if multiple_max[1] == 'd':
                print('eerste en tweede hetzelfde')
                # two bins are taken together, random phase ratio doubles 9.6%
                highest = (int(multiple_max[0]) + int(multiple_max[0]))
                phase_ratio_t = highest/total

                twobin= bindf.Reads[bindf.Reads == int(multiple_max[0])].index.tolist()
                HB_two = bindf.copy()
                HB_two['HeadBin'] = 'N'
                for x in twobin:
                    HB_two = HB_two.set_value(x, 'HeadBin', 'Y')
                print(phase_ratio_t)

                # random chanche/bin 4.8%, highest distinct reads is taken
                highest_d = bindf['DistinctReads'].max() #highest distinct
                onebin= bindf.DistinctReads[bindf.DistinctReads == highest_d].index.tolist()
                HB_one = bindf.copy()
                HB_one['HeadBin'] = 'N'
                for x in onebin:
                    HB_one = HB_one.set_value(x, 'HeadBin', 'Y')
                print(HB_one)
                phase_ratio_o = highest_d/total
                print(phase_ratio_o)

            elif multiple_max[2] == 'd':
                print('tweede en derde hetzelfde')
                # random chanche/bin 4.8%, highest reads is taken
                highest = bindf['Reads'].max()
                phase_ratio_o = highest/total
                onebin= bindf.Reads[bindf.Reads == highest].index.tolist()
                HB_one = bindf.copy()
                HB_one['HeadBin'] = 'N'
                for x in onebin:
                    HB_one = HB_one.set_value(x, 'HeadBin', 'Y')
                print(phase_ratio_o)

                # two highest bins are taken
                dropmax = bindf[bindf.Reads != multiple_max[0]]          #to take from the second bin, the highest distinct, 
                                                                        #the highest bin is dropped
                highest_d = dropmax['DistinctReads'].max() #highest distinct
                #starti = np.where(dropmax['DistinctReads'] == highest_d) #index of highest distinct
                #print(starti)
                #starti = str(starti)
                #print(starti)
                #starti = starti.replace('(array([', '')
                #print(starti)
                #starti = starti.replace(']),)', '')
                #print(starti)
                #starti = int(starti) + 1
                #print(starti)
                twobin= bindf.Reads[bindf.Reads == int(multiple_max[0])].index.tolist()
                twobin.append(int(starti))
                HB_two = bindf.copy()
                HB_two['HeadBin'] = 'N'
                for x in twobin:
                    HB_two = HB_two.set_value(x, 'HeadBin', 'Y')
                print(HB_two)
                phase_ratio_t = (int(multiple_max[0])+int(multiple_max[1]))/total
                print(phase_ratio_t)

            else:
                print('drie verschillenden')
                highest = bindf['Reads'].max()
                phase_ratio_o = highest/total
                onebin= bindf.Reads[bindf.Reads == highest].index.tolist()
                HB_one = bindf.copy()
                HB_one['HeadBin'] = 'N'
                for x in onebin:
                    HB_one = HB_one.set_value(x, 'HeadBin', 'Y')
                print(phase_ratio_o)

                # two highest bins are taken
                dropmax = bindf[bindf.Reads != multiple_max[0]]          #to take from the second bin, the highest distinct, 
                                                                        #the highest bin is dropped
                highest_d = dropmax['DistinctReads'].max() #highest distinct
                #starti = np.where(dropmax['DistinctReads'] == highest_d) #index of highest distinct
                #starti = str(starti)
                #starti = starti.replace('(array([', '')
                #starti = starti.replace(']),)', '')
                #starti = int(starti) + 1
                twobin= bindf.Reads[bindf.Reads == int(multiple_max[0])].index.tolist()
                twobin.append(int(starti))
                HB_two = bindf.copy()
                HB_two['HeadBin'] = 'N'
                for x in twobin:
                    HB_two = HB_two.set_value(x, 'HeadBin', 'Y')
                phase_ratio_t = (int(multiple_max[0])+int(multiple_max[1]))/total
                print(phase_ratio_t)

        elif length == 1:
            print('not unique high bin')

        total = bindf['Reads'].sum()
        total2 = bindf['DistinctReads'].sum()
        HB_two = HB_two.set_index('BIN').T
        HB_one = HB_one.set_index('BIN').T
        HB_one.replace('', 'NA')
        HB_two['Sum'] = [total, total2, 'NA']
        HB_one['Sum'] = [total, total2, 'NA']
        HB_two['PhaseRatio'] = [phase_ratio_t, 'NA', 'NA']
        HB_one['PhaseRatio'] = [phase_ratio_o, 'NA', 'NA']
        #print(HB_one)
        #print(HB_two)
        #print(directory)
        out2name = str(directory + '/' + fname + '_1Bin.txt')
        out3name = str(directory + '/' + fname + '_2Bin.txt')
        #print(out2name)
        #print(out3name)
        HB_one.to_csv(out2name, index=False, sep="\t")
        HB_two.to_csv(out3name, index=False, sep="\t")
        #return (df, HB_one)

def score_vs_ratio(self, bin_path):
    path1 = self.shs_path
    acc = os.path.basename(os.path.dirname(path1))
    print(acc)
    pr_dict = {}
    
    for result in glob.iglob(bin_path):
        name = result
        name = os.path.basename(result)
        name = name.replace('Bin.txt', '')
        df = pd.read_csv(result, sep = '\t')
        pr = df.amount[22]
        pr_dict[name] = pr
    df1 = pd.read_csv(path1, sep = '\t')
    cluster_pos = df1.columns.get_loc("Name")
    
    PR_list = []
    for cluster in df1['Name']:
        PR = pr_dict.get(cluster)
        PR_list.append(PR)
    insert_pos = (df1.columns.get_loc("PhaseScore"))+1
    PRnum = str('PhasedRatio' + str(m))
    df1.insert(insert_pos, PRnum, PR_list)
    out2name = str('/Users/michelle/'+ acc + 'wPR.txt')   #aanpassen
    df1.to_csv(out2name, index=False, sep="\t")
