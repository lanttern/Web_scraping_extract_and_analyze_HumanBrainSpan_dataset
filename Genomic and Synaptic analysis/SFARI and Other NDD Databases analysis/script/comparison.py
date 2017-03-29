# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 09:07:11 2014

@author: zhihuixie
"""

import pandas as pd
from Tkinter import *
from tkFileDialog import askopenfilename
import scipy as sp
import scipy.stats

class Comparison():
    """
    """
    def __init__(self, list1, list2):
        """
        """
        self.list1 = list1
        self.list2 = list2
        
    def compare(self):
        """
        """
        com = [item for item in self.list1 if item in self.list2]
        return len(com), com
        
    def paired_number(self, total_num):
        """
        """
        com_num, _ = self.compare()
        return [com_num, total_num - com_num]
        
    def fisher_exact_test(self, pairs):
        """
        """
        percentage_1 = float(pairs[0][0])*100/(pairs[0][0]+ pairs[0][1])
        percentage_2 = float(pairs[1][0])*100/(pairs[1][0] + pairs[1][1])
        print "pairs:", pairs
        print "Percentage of NDDs: %.2f %.2f"%(percentage_1, percentage_2) 
        print "Enrichment factor: %.2f" %(percentage_1/percentage_2)
        print "Result of fisher exact test: %.2E"%(sp.stats.fisher_exact(pairs)[1]*28)
        
fname = "unsigned"
dfs = []
def openFile():
    """
    """
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
        
    list1 = dfs[0].Symbol.values
    list2 = dfs[1].Symbol.values
    list3 = dfs[2].Symbol.values
    
    #compare enrichment of NDDs associated MET Co-IP candidates at genomic levels
    #list1 is 72 candaidates
    #list2 is SFARI candidates of category 1 and category 2
    compare  = Comparison(list1, list2)
    #NDDs associated candidates in MET Co-IP
    com_num, com = compare.compare()
    print "NDDs associated candidates in MET Co-IP in genome:"
    print com_num, com, "\n"
    #pair of MET Co-IP category
    pair1 = compare.paired_number(len(list1))
    #pair of Genome category
    pair2 = [len(list2), 30000 - len(list2)]
    #fisher exact test
    print "Comparison of NDDs candidates in MET Co-IP at genomic levels:"
    compare.fisher_exact_test([pair1, pair2])
    print "\n"

    #compare enrichment of NDDs associated MET Co-IP candidates at synaptic levels
    list4 = com
    #list1 is 72 candaidates
    #list2 is SFARI candidates of category 1 and category 2
    #list3 is synaptosome proteins
    #list4 is NDDs associated candidates in MET Co-IP
    #candidates in MET Co-IP and synaptosome
    compare0  = Comparison(list3, list1)
    com_num_syn_MET, com_syn_MET = compare0.compare()
    print "Candidates in MET Co-IP and synaptosome:"
    print com_num_syn_MET, com_syn_MET, "\n"
    #NDDs associated candidates in MET Co-IP and synaptosome
    compare1  = Comparison(list3, list4)
    com_num_syn_MET_ndd, com_syn_MET_ndd = compare1.compare()
    print "NDDs associated candidates in MET Co-IP in synaptosome:"
    print com_num_syn_MET_ndd, com_syn_MET_ndd, "\n"
    #NDDs associated candidates in synaptosome
    compare2  = Comparison(list3, list2)
    com_num_syn_ndd, com_syn_ndd = compare2.compare()
    print "NDDs associated candidates in synaptosome:"
    print com_num_syn_ndd, com_syn_ndd, "\n"
    #pair of MET Co-IP category
    pair3 = compare1.paired_number(com_num_syn_MET)
    #pair of Genome category
    pair4 = [com_num_syn_ndd, len(list3) - com_num_syn_ndd]
    #fisher exact test
    print "Comparison of NDDs candidates in MET Co-IP at genomic levels:"
    compare.fisher_exact_test([pair3, pair4])
    print "\n"