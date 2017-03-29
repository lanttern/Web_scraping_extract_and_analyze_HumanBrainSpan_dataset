# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 20:42:43 2014

@author: zhihuixie
"""

import pandas as pd
import re


### the following lines are for Syn_DB counting
df0 = open("SynDB_mouse_database.csv")
mwords = ""
for row in df0:
    mwords += row

mstart_index = [m.start() for m in re.finditer("Name:", mwords)]
mend_index = [m.start() for m in re.finditer(" Description", mwords)]
ms1 = [m.start() for m in re.finditer("_MOUSE", mwords)]
mgene_symbols = []
muniport_ids = []
for i in range (len(mstart_index)):
    mgene_symbol = mwords[mstart_index[i] + 5:mend_index[i]]
    mgene_symbol = mgene_symbol.upper()
    if len(mgene_symbol) != 0:
        mgene_symbols.append(mgene_symbol)
    else:
        #print start_index[i], end_index[i], "test"
        #print s1[i-1], e1[i]
        muniport_id = mwords[ms1[i] - 6: ms1[i]]
        muniport_ids.append(muniport_id)
df_mouse = pd.DataFrame({"Gene symbols": mgene_symbols})
df_mouse.to_csv("Mouse SynDB dataset.csv", index = False)
df = open("SynDB_human_database.csv")
words = ""
for row in df:
    words += row

start_index = [m.start() for m in re.finditer("Name:", words)]
end_index = [m.start() for m in re.finditer(" Description", words)]
s1 = [m.start() for m in re.finditer("_HUMAN", words)]
gene_symbols = []
uniport_ids = []
for i in range (len(start_index)):
    gene_symbol = words[start_index[i] + 5:end_index[i]]
    if len(gene_symbol) != 0:
        gene_symbols.append(gene_symbol)
    else:
        #print start_index[i], end_index[i], "test"
        #print s1[i-1], e1[i]
        uniport_id = words[s1[i] - 6: s1[i]]
        uniport_ids.append(uniport_id)
df_human = pd.DataFrame({"Gene symbols": gene_symbols})
df_human.to_csv("Human SynDB dataset.csv", index = False)        

"""
# the following lines are for hPSD proteomics database counting
df = pd.read_csv("hPSD_proteomics.csv")
gene_symbols = df.Symbol.values
print gene_symbols
#print gene_symbols, len(uniport_ids), len(gene_symbols)
#print file1, file1[0]
#class Count (files):
candidates = []
df1 = open("Candidates_72.csv")
for line in df1:
    line = line.split()
    #line = line.replace("\n","")
    candidates.append(line[0])

synaptic_candidates = [item for item in candidates if item in gene_symbols]

c1 =[]
df2 = open("Candidates_15.csv")
for line1 in df2:
    line1 = line1.split()
    #line = line.replace("\n","")
    c1.append(line1[0])
print c1

ndd_candidates = [item for item in c1 if item in gene_symbols]

df3 = open("sfari_ndd.csv")

for line2 in df3:
    line2 = line2.split()
    #line = line.replace("\n","")


ndd_synaptic_genes = [item for item in gene_symbols if item in line2]

df4 = open("NDDCandidatesTotalNetwork-11.csv")
c2 = []
for line3 in df4:
    line3 = line3.split()

t_ndds = [item for item in line3 if item in gene_symbols]
print t_ndds
df5 = open("NDDCandidatesSubNetwork-12.csv")

for line4 in df5:
    line4 = line4.split()

s_ndds = [item for item in line4 if item in gene_symbols]

df6 = pd.read_csv("FullBuild_nodes.csv")
total_nodes = df6.Symbol.values
print len([item for item in total_nodes if item in gene_symbols])

df7 = pd.read_csv("Subnetwork_nodes.csv")
sub_nodes = df7.Symbol.values
print len([item for item in sub_nodes if item in gene_symbols])

print len(ndd_synaptic_genes)


print "candidates in full network associated with ndd: %d"%len(t_ndds)
print "candidates in sub network associated with ndd: %d"%len (gene_symbols)
print "synaptic candidates in MET pulldown associated with ndd: %d"%len(ndd_candidates)
print "synaptic candidates in MET pulldown: %d"%len(synaptic_candidates)
print "synaptic candidates associated with ndd: %d"%len(ndd_synaptic_genes)
print len(s_ndds)
#print float(len([a for a in mgene_symbols if a in gene_symbols]))/len(mgene_symbols)
#print [a for a in synaptic_candidates if a in mgene_symbols]
"""