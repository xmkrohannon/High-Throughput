#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 18:39:23 2019

@author: xander
"""

import matplotlib.pyplot as plt, os, pandas as pd
os.path.join(os.path.dirname(__file__))

#Open config file to create list of all genes, in order, contained in all solutions
gene_List = []
with open ('TF-Gene.cfg') as file:
    text = file.readlines()
    text = text[31:1905]
    for line in text:
        line = line.strip().split('\t')
        gene_List.append(line[1])
file.close()

#Read all solution files and gather values for all genes in the correct order
model_Data = []
exp_Names = []
file_Names = ['TF-Gene_solution_4.dat','TF-Gene_solution_8.dat','TF-Gene_solution_10.dat','Omics_solution_10.dat']
for name in file_Names:
    with open (name) as file:
        text = file.readlines()
        for i in range (0,len(text)):
            line = text[i]
            line = line.strip().split('\t')
            line = line[2:]
            count = 0
            for j in range(0,len(line),len(gene_List)):
                steady_State = list(map(float,line[j:j+len(gene_List)]))
                model_Data.append(steady_State)
                count += 1
                title = name.split('.')[0] + '_' + str(i + 1) + '_' + str(count)
                exp_Names.append(title)

#Create data frame associating experiments with gene states
data_Frame = pd.DataFrame(model_Data, columns=gene_List, index=exp_Names).transpose()

#Perform k-means biclustering on the data frame
models = []
scores = []
clusters = []
from sklearn import cluster
for i in range (2,int(data_Frame.shape[1]**(0.5))):
    clusters.append(i)
    model = cluster.bicluster.KMeans(n_clusters=i)
    model.fit(data_Frame)
    models.append(model)
    scores.append(model.inertia_)

k_Vals = list(range(2,int(data_Frame.shape[1]**(0.5))))

#Create plot for elbow method for determining ideal bicluster
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(k_Vals, scores, 'bx-')
plt.xlabel('k')
plt.ylabel('Sum_of_squared_distances')
plt.title('Elbow Method For Optimal k')
plt.show()
fig.savefig('K-means_Scores.png')

#Export data frame to csv format for use in Morpheus
data_Frame.to_csv('Steady_State_Analysis.csv', sep=',')