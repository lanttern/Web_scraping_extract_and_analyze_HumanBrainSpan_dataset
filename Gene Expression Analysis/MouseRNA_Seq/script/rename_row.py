# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 14:19:50 2014

@author: zhihuixie
"""
import pandas as pd
df = pd.read_csv("mouse_rna_seq.csv", index_col = 0)
a = []
for index, row in df.iterrows():
    #print index
    index1 = index[19:]
    a.append(index1)
df1 = pd.DataFrame(a)
df1.to_csv("test1.csv")