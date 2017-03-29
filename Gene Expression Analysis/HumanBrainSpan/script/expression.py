# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:39:34 2014

@author: zhihuixie
"""

import pandas as pd
from Tkinter import *
from tkFileDialog import askopenfilename

class Expression_data():
    def __init__(self, gene_list, data):
        """
        """
        self.gene_list = gene_list
        self.data = data
        
    def write_to_file(self, file_dir):
        """
        """
        extracted_data = self.data[self.data["Symbol"].isin(self.gene_list)]
        extracted_data.to_csv(file_dir, index = False)

fname = "unsigned"
dfs = []
def openFile():
    """
    """
    global fname
    fname = askopenfilename(parent = root)
    root.destroy()

if __name__ == '__main__':
    for i in range (2):
        root = Tk()
        root.update()
        Button(root, text='File Open', command = openFile).pack(fill=X)
        mainloop()
        dfs.append(pd.read_csv(fname))
        
#open MET 72 candidates first
df_MET_candidates = dfs[0]
gene_list = df_MET_candidates.Symbol.values
print gene_list
#open RNA seq data
df_rna_seq = dfs[1]
#save extracted data to file
file_dir = "../data/expression data/clean data/gene expression of MET-interacting candidates-S1C.csv"
extracted_data = Expression_data(gene_list, df_rna_seq)
extracted_data.write_to_file(file_dir)
