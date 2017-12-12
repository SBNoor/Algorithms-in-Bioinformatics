#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:23:08 2017

@author: Noor
"""

import re
import numpy as np
from Bio import Phylo
import sys

def extracting_distance_matrix(filename):
    names = []
    lst_of_dist = []
    
    with open(filename) as f:
        line = f.readline()
        line = f.readline()
        
        while line:
            x = re.split('\s+',line.strip())
            names.append(x[0])
            lst = []
            for i in range(1,len(x)):
                lst.append(float(x[i]))
            lst_of_dist.append(lst)
            line = f.readline()
            
    distance_matrix = np.zeros(shape=(len(lst_of_dist),len(lst_of_dist[0])))
    for i in range(len(lst_of_dist)):
        for j in range(len(lst_of_dist[i])):
            distance_matrix[i][j] = lst_of_dist[i][j]
    return distance_matrix,names


def nj(distance_matrix,names,count):
    #r = []
    r = [0] * len(distance_matrix)
    while (len(distance_matrix)) > 3:
        r = average_dist_to_all_tips(distance_matrix,names,r)
        row,column = finding_lowest_value(distance_matrix,r)
        count = branch_length(distance_matrix,row,column,r,count)
        distance_matrix = update_dist_matrix(distance_matrix, row, column)
        
    return termination(count)

def average_dist_to_all_tips(distance_matrix,names,r):
    for i in range(len(distance_matrix)):
        #r.append(0)
        r[i] = 0
        for j in range(len(distance_matrix[0])):
            r[i] += distance_matrix[i][j]
        r[i] = r[i]/(len(distance_matrix)-2)     
    return r

def finding_lowest_value(distance_matrix,r):
    closest_pair_value = distance_matrix[1][0] - (r[1] + r[0])
    row = 1
    column = 0
    for i in range(1,len(distance_matrix)):
        for j in range(len(distance_matrix)):
            if ((distance_matrix[i][j] - (r[i] + r[j]))<closest_pair_value) and (i != j):
                closest_pair_value = distance_matrix[i][j] - (r[i] + r[j])
                row = i
                column = j
                
    return(row,column)
    
def update_dist_matrix(distance_matrix, row, col):
    for k in range(0, len(distance_matrix)):
        if k != row and k != col:
            distance_matrix[col][k] = (distance_matrix[row][k] + distance_matrix[col][k] - distance_matrix[row][col]) / 2
            distance_matrix[k][col] = distance_matrix[col][k]               

    distance_matrix = np.delete(distance_matrix,(row),axis = 1)
    distance_matrix = np.delete(distance_matrix,(row),axis = 0)
    return distance_matrix
        
def branch_length(distance_matrix,row,col,r,count):
    v1 = v2 =  0
    
    v1 = (distance_matrix[row][col] + (r[row] - r[col]))/2
    v2 = (distance_matrix[row][col] + (r[col] - r[row]))/2
    
    #print("v1 : ",v1," v2 : ",v2)
    
    clade1 = clades[row];
    clade2 = clades[col];
    count += 1;
    inner_clade = Phylo.BaseTree.Clade(None, "Node " + str(count));
    inner_clade.clades.append(clade1);
    inner_clade.clades.append(clade2);
    
    # assign branch length
    clade1.branch_length = v1
    clade2.branch_length = v2
                            
    # update node list
    clades[col] = inner_clade;
    del clades[row];
    
    return count

def termination(count):
    count += 1;
    clade_k = Phylo.BaseTree.Clade(None, "Node " + str(count));
    clade1 = clades[0];
    clade2 = clades[1];
    clade3 = clades[2];
    
    clade_k.clades.append(clade1);
    clade_k.clades.append(clade2);
    clade_k.clades.append(clade3);
    
    # assigning branch length
    clade1.branch_length = (distance_matrix[0, 1] + distance_matrix[0, 2] - distance_matrix[1, 2]) / 2.0;
    clade2.branch_length = (distance_matrix[0, 1] + distance_matrix[1, 2] - distance_matrix[0, 2]) / 2.0;
    clade3.branch_length = (distance_matrix[0, 2] + distance_matrix[1, 2] - distance_matrix[0, 1]) / 2.0;

    # Constructing the tree
    tree = Phylo.BaseTree.Tree(clade_k, rooted=False);
    return tree

if __name__=="__main__":

    filename = sys.argv[1]
    distance_matrix,names = extracting_distance_matrix(filename)
    count = 0
    clades = [Phylo.BaseTree.Clade(None, name) for name in names];
    tree = nj(distance_matrix,names,count)
    out_name = str(filename[:-4]).split("/")[1]
    Phylo.write(tree, out_name, 'newick')
    
    # tree = Phylo.read('output-tree.new', 'newick')
    # Phylo.draw(tree)

    
    
    
