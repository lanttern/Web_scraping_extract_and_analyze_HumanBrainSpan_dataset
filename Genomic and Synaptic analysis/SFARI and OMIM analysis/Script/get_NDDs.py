# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:14:08 2015

@author: zhihuixie
"""

from Tkinter import *
from tkFileDialog import askopenfilename
import pandas as pd

fname = "unsigned"
dfs = []
def openFile():
    global fname
    fname = askopenfilename(parent = root)
    root.destroy()

if __name__ == '__main__':
    for i in range (3):
        root = Tk()
        root.update()
        Button(root, text='File Open', command = openFile).pack(fill=X)
        mainloop()
        dfs.append(pd.read_csv(fname))
    df_sfari_score = dfs[0]
    df_MET_candidates_72 = dfs[1]
    df_MET_candidates_full_network = dfs[2]
    sfari_candidates = df_sfari_score["Gene Symbol"].values
    MET_candidates_72 = df_MET_candidates_72["Gene symbols"].values
    MET_candidates_full_network = df_MET_candidates_full_network["Candidates"].values
    print len(sfari_candidates)
    MET_candidates_72_NDDs = [candidate for candidate in MET_candidates_72 \
                             if candidate in sfari_candidates]
                                 
    MET_candidates_full_network_NDDs = [candidate for candidate in \
                                        MET_candidates_full_network \
                                        if candidate in sfari_candidates]
    print len(MET_candidates_72_NDDs), MET_candidates_72_NDDs
    print len(MET_candidates_full_network_NDDs), MET_candidates_full_network_NDDs
        