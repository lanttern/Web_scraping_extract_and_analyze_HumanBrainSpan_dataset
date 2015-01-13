# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:39:34 2014

@author: zhihuixie
"""

import pandas as pd
import csv
class Expression_data():
    def __init__(self, file1, file2, structure):
        self.file1 = file1
        self.file2 = file2
        self.structure = structure
        self.data = open(self.file1)
        self.df = pd.read_csv(self.file2)
    def extract_expression_data(self):
        print self.data
        for candidate in self.data:
            candidate = candidate.replace("\r", " ")
            candidate = candidate.split()
        candidate = ["age"] + candidate
        print candidate
        return candidate
    def write_all_to_file(self):
        self.genes = self.extract_expression_data()
        temp = self.df[self.df["gene_symbol"].isin(self.genes)]
        temp = temp.drop("Unnamed: 0", axis = 1)
        temp.to_csv(self.file1[:-4] + "_all_structures_expression.csv")
    def write_structure_to_file(self):
        column_values = self.df.columns.values.tolist()
        structures = ["gene_symbol"]
        for item in self.structure:
            for value in column_values:
                if item in value:
                    structures.append(value)
        new_data = self.df.ix[:, structures]
        new_data = new_data[new_data["gene_symbol"].isin(self.genes)]
        new_data.to_csv(self.file1[: -4] + "_" + self.structure[0] + "_to_" + self.structure[-1] + "_expression.csv")
file1 = "MET.csv"
file2 = "expression_matrix_add_Rows_and_Columns_full_version.csv"
STC_data = Expression_data(file1,file2, ["STC", "ITC", "TCx", "V1C"])
STC_data.write_all_to_file()
STC_data.write_structure_to_file()
