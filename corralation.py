# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 16:19:33 2014

@author: zhihuixie
"""

import pandas as pd

def extract_expression_data(file1):
    read_file = open(file1)
    for candidate in read_file:
        candidate = candidate.replace("\r", " ")
        candidate = candidate.split()
    return candidate

def creat_correlation_table(file1):
    genes = extract_expression_data(file1)
    #print genes
    output = pd.read_csv("1.csv")
    output = output[output["gene-symbol"].isin (genes)]
    for i in range (2, 28):
        temp_file = str(i) + ".csv"
        df = pd.read_csv(temp_file)
        d = df[df["gene-symbol"].isin (genes)]
        output = pd.merge(output, d, how='outer')
    """       
    for f in files:
        print f
        output.merge(f)
    """
   # print output

    output.to_csv("result.csv")
    
creat_correlation_table("NDDCandidatesSubNetwork-12.csv")