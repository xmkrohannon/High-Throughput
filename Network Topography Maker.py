#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 13:20:41 2019

@author: xander
"""

import os
os.path.join(os.path.dirname(__file__))

def TopoMaker(input_File):
    acting_Genes = []
    target_Genes = []
    gene_Action = []
    with open (input_File) as file:
        text = file.readlines()
        for i in range (35,len(text)):
            line = text[i].strip().split('\t')
            acting_Gene = line[0]
            if (acting_Gene == 'H-NS'):
                acting_Gene = 'hns'
            elif (acting_Gene == 'IHF'):
                acting_Gene = 'ihfa-ihfb'
            elif (acting_Gene == 'HU'):
                acting_Gene = 'hupA-hupB'
            elif (acting_Gene == 'HipAB'):
                acting_Gene = 'hipA-hipB'
            elif (acting_Gene == 'RcsAB'):
                acting_Gene = 'rcsA-rcsB'
            if ('-' in acting_Gene):
                acting_Gene = acting_Gene.split('-')
                gene_One = acting_Gene[0]
                if (len(gene_One) == 3):
                    gene_One = gene_One.lower()
                else:
                    gene_One = gene_One[:1].lower() + gene_One[1:]
                gene_Two = acting_Gene[1]
                if (len(gene_Two) == 3):
                    gene_Two = gene_Two.lower()
                else:
                    gene_Two = gene_Two[:1].lower() + gene_Two[1:]
                target_Gene = line[1]
                action = line[2]
                if (len(action) == 2):
                    acting_Genes.append(gene_One)
                    target_Genes.append(target_Gene)
                    gene_Action.append('1')
                    acting_Genes.append(gene_One)
                    target_Genes.append(target_Gene)
                    gene_Action.append('2')
                    acting_Genes.append(gene_Two)
                    target_Genes.append(target_Gene)
                    gene_Action.append('1')
                    acting_Genes.append(gene_Two)
                    target_Genes.append(target_Gene)
                    gene_Action.append('2')
                else:
                    if (action == '+'):
                        acting_Genes.append(gene_One)
                        target_Genes.append(target_Gene)
                        gene_Action.append('1')
                        acting_Genes.append(gene_Two)
                        target_Genes.append(target_Gene)
                        gene_Action.append('1')
                    elif (action == '-'):
                        acting_Genes.append(gene_One)
                        target_Genes.append(target_Gene)
                        gene_Action.append('2')
                        acting_Genes.append(gene_Two)
                        target_Genes.append(target_Gene)
                        gene_Action.append('2')
            else:
                if (len(acting_Gene) == 3):
                    acting_Gene = acting_Gene.lower()
                else:
                    acting_Gene = acting_Gene[:1].lower() + acting_Gene[1:]
                target_Gene = line[1]
                action = line[2]
                if (len(action) == 2):
                    acting_Genes.append(acting_Gene)
                    target_Genes.append(target_Gene)
                    gene_Action.append('1')
                    acting_Genes.append(acting_Gene)
                    target_Genes.append(target_Gene)
                    gene_Action.append('2')
                else:
                    if (action == '+'):
                        acting_Genes.append(acting_Gene)
                        target_Genes.append(target_Gene)
                        gene_Action.append('1')
                    elif (action == '-'):
                        acting_Genes.append(acting_Gene)
                        target_Genes.append(target_Gene)
                        gene_Action.append('2')
    file.close()
    topo_Data = [acting_Genes,target_Genes,gene_Action]
    
    #Used for creating the TF-Gene data file
    '''with open ('TF-Gene.topo','w') as file:
        file.write('Source' + '\t' + 'Target' + '\t' + 'Type' + '\n')
        for i in range (0,len(topo_Data[0])):
            file.write(topo_Data[0][i] + '\t' + topo_Data[1][i] + '\t' + topo_Data[2][i] + '\n')    
    return topo_Data
    file.close()'''
    
    '''#Used for adding TF-TF to the combined data file
    with open ('Omics.topo','a') as file:
        for i in range (0,len(topo_Data[0])):
            file.write(topo_Data[0][i] + '\t' + topo_Data[1][i] + '\t' + topo_Data[2][i] + '\n')
    file.close()'''

#TopoMaker('TF_Gene Interaction Network.txt')
#TopoMaker('TF_TF Interaction Network.txt')