# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 19:05:23 2014

@author: zhihuixie
"""

from bs4 import BeautifulSoup
import urllib2
import csv
file1 = urllib2.urlopen("http://api.brain-map.org/api/v2/data/query.xml?criteria=service::dev_human_correlation[set$eqrna_seq_genes][probes$eq1090294][structures$eq'CBC'][num_rows$eq500]")
file2 = open("test.csv", "wb")
fields = ["id","name","gene-id","gene-symbol","gene-name", "entrez-id","chromosome", "start", "end", "r"]
wr = csv.DictWriter(file2, fieldnames = fields)
wr.writeheader()
soup = BeautifulSoup(file1)
soup_object = soup.find_all("object")
print soup_object[-1].find("gene-symbol").text
data = []
for object1 in soup_object:
    #wr.writerow(object1.text)
    row = {}
    for item in fields:
        if type(object1.find(item)) == type(None):
            #print object1.find(item), item
            continue
        row[item] = str(object1.find(item).text)
    wr.writerow(row)
    #row.append(object1.text)
    #print row
    #break