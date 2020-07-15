#!/.usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:04:58 2019

@author: xander
"""

import copy, os, pandas as pd
os.path.join(os.path.dirname(__file__))

#Open biclustered file from Morpheus
with open ('Cluster_Steady_State.gct') as file:
    text = file.readlines()
    cols = text[2].strip().split('\t')
    cols.pop(0)
    rows = []
    data = []
    for i in range(3,len(text)):
        line = text[i]
        line = line.strip().split('\t')
        data.append(line)
        rows.append(data[-1].pop(0))
file.close()

#Create data frame from file
data_Frame = pd.DataFrame(data, columns = cols, index = rows)
data_Frame = data_Frame.replace('na', 0)
data_Frame = data_Frame.astype('float')

#Create dictionary of experiments contained in each cluster
sim_Clusters = {}
for i in range (1,16):
    for j in range (1, len(cols)):
        if (data_Frame[cols[j]][0] == i):
            if (i not in sim_Clusters.keys()):
                sim_Clusters[i] = []
            sim_Clusters[i].append(cols[j])

#Check if TF-Gene solutions and TF-TF-Gene solutions were clustered together
#See what modes of solutions were clustered together
co_Clusters = []
cls_Clusters = {}
for i in sim_Clusters.keys():
    hits = []
    tgen_Chk = False
    omic_Chk = False
    for j in sim_Clusters[i]:
        if ('TF-Gene' in j):
            tgen_Chk = True
        if ('Omics' in j):
            omic_Chk = True
        line = j.split('_')
        hits.append(line[2])
    if (tgen_Chk and omic_Chk):
        co_Clusters.append(i)
    cls_Clusters[i] = list(set(hits))

#Create a copy of the original data frame, remove the cluster IDs
#and calculate the row, column, and overall means for Cheng and Church method
rev_Frame = copy.deepcopy(data_Frame)
row_Label = rev_Frame[cols[0]]
col_Label = rev_Frame.loc[rows[0]]
rev_Frame = rev_Frame.drop(cols[0])
rev_Frame = rev_Frame.drop(columns = rows[0])
total_Mean = rev_Frame.stack().mean()
curr_Mean = rev_Frame.stack().mean()
col_Mean = rev_Frame.mean(axis = 0)
row_Mean = rev_Frame.mean(axis = 1)

#Perform Cheng and Church method for simplifying bicluster
alpha = 1
beta = 0.2
for i in range (0,rev_Frame.shape[0]):
    if (row_Mean[i] > alpha*total_Mean):
        rev_Frame = rev_Frame.drop(rows[i + 1])
        row_Label = row_Label.drop(rows[i + 1])
        curr_Mean = rev_Frame.stack().mean()
    elif ((curr_Mean - rev_Frame.drop(rows[i + 1]).stack().mean()) > beta) and (curr_Mean > beta):
        rev_Frame = rev_Frame.drop(rows[i + 1])
        row_Label = row_Label.drop(rows[i + 1])
        curr_Mean = rev_Frame.stack().mean()
for j in range (0, rev_Frame.shape[1]):
    if (col_Mean[j] > alpha*total_Mean):
        rev_Frame = rev_Frame.drop(columns = cols[j + 1])
        col_Label = col_Label.drop(cols[j + 1])
        curr_Mean = rev_Frame.stack().mean()
    elif ((curr_Mean - rev_Frame.drop(columns = cols[j + 1]).stack().mean()) > beta) and (curr_Mean > beta):
        rev_Frame = rev_Frame.drop(columns = cols[j + 1])
        col_Label = col_Label.drop(cols[j + 1])
        curr_Mean = rev_Frame.stack().mean()

#Reorder and export revised bicluster for Morpheus
rev_Frame.insert(0, cols[0], row_Label)
rev_Rows = rev_Frame.index
rev_Rows = rev_Rows.insert(0, rows[0])
rev_Frame = rev_Frame.append(col_Label)
rev_Frame = rev_Frame.loc[rev_Rows]
rev_Frame.to_csv('Revised_Steady_State_Analysis.csv', sep=',')

#Create dictionary of genes contained in each cluster
gene_Clustsers = {}
for i in range (1,16):
    for j in range (1, len(rev_Rows)):
        if (rev_Frame[cols[0]][j] == i):
            if (i not in gene_Clustsers.keys()):
                gene_Clustsers[i] = []
            gene_Clustsers[i].append(rev_Rows[j])

#Create a file for gene ontologies using Panther
for i in gene_Clustsers.keys():
    with open ('Gene_Cluster' + str(i) + '.txt', 'w') as file:
        for j in range(0,len(gene_Clustsers[i])):
            if (j < (len(gene_Clustsers[i]) - 1)):
                file.write(gene_Clustsers[i][j] + '\n')
            else:
                file.write(gene_Clustsers[i][j])
    file.close()