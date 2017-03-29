# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:39:34 2014

@author: zhihuixie
"""

import pandas as pd
class Expression_data():
    def __init__(self, file1, file2):
        self.file1 = file1
        self.file2 = file2
        self.data = open(self.file1)
        self.df = pd.read_csv(self.file2)
    def extract_expression_data(self):
        for candidate in self.data:
            candidate = candidate.replace("\r", " ")
            candidate = candidate.split()
        candidate = ["Geneid"] + candidate
        print candidate
        return candidate
    def write_to_file(self):
        genes = self.extract_expression_data()
        temp = self.df[self.df["Geneid"].isin(genes)]
        temp.to_csv(self.file1[:-4] + "_output.csv")
    def write_all_to_file(self):
        self.df.to_csv(self.file2 + "_output.csv")
        list1 = ["STC." + str(num) for num in range(1,36)]
        list1 = ["STC"] + list1
        new_data = self.df.ix[:, list1]
        new_data.to_csv(self.file2[: -4] + "_output.csv")
#file1 = "Candidates_15.csv"
#file2 = "expression_matrix_addRowand_Column.csv"
#STC_data = Expression_data(file1,file2)
#STC_data.write_all_to_file()

#STC.35,ITC.33
#file1 = "Candidates_15.csv"
#file2 = "BrainSpan_expression_matrix.csv"
file3 = "Candidates_72.csv"
file5 = "mouse_rna_seq_rename.csv"
#file4 = "NDDCandidatesTotalNetwork-11.csv"
#NDDs_in_MET_pulldown = Expression_data(file1, file2)
#NDDs_in_SubNetwork = Expression_data(file3, file2)
#NDDs_in_TotalNetwork = Expression_data(file4, file2)
#NDDs_in_MET_pulldown.write_to_file()
#NDDs_in_SubNetwork.write_to_file()
#NDDs_in_TotalNetwork.write_to_file()
NDDs_in_SubNetwork_mouse = Expression_data(file3, file5)
NDDs_in_SubNetwork_mouse.write_to_file()
