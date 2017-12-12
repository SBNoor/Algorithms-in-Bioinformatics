#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import Phylo
import numpy as np
import sys



def radixSortFixedString(array):
    fixedLength = len(array[0])
    oa = ord('0'); 
    oz = ord('1'); 
    n = oz - oa + 1; 
    buckets = [[] for i in range(0, n)] 
    for position in reversed(range(0, fixedLength)):
        for string in array:
            buckets[ord(string[position]) - oa].append(string) 
        del array[:]
        for bucket in buckets: # Reassemble array in new order
            array.extend(bucket)
            del bucket[:]
    return array


def create_bit_vector(lst, size):
    
    bit_vector = np.ones(size)
    for j in range(len(lst)):
        number = lst[j].name.split('@')[0]
        index = int(number)
        bit_vector[index -1] = 0
        
    return [int(i) for i in bit_vector]
    


def tabulate_names_one(tree):
    names = {}
    iter = 0
    for idx, clade in enumerate(tree.find_clades(order = "postorder")):
        if clade.name:
            clade.name = '%d@%s' % (iter, clade.name)
            iter +=1
        else:
            clade.name = 'innernode'
        names[clade.name] = clade
    return names

def tabulate_names_two(tree,lst):
    names = {}
    for idx, clade in enumerate(tree.find_clades(order = "postorder")):
        if clade.name:
            for i in range(len(lst)):
                name = lst[i].name.split('@')[1]
                if name == clade.name:
                    clade.name = lst[i].name
                    break

        else:
            clade.name = 'innernode'
        names[clade.name] = clade
    return names

def get_look_up_table(tree):
    for clade in tree.get_nonterminals(order = "postorder"):
        lst = clade.get_terminals(order = "postorder")
    return lst

def min_node(lst):
    
    numbers = []

    for i in range(len(lst)):
        name = lst[i].name.split('@')[0]

        
        numbers.append(int(name))
    min = numbers[0]
    for i in numbers:
        if min >= i:
            min = i
    return min
    
    
def max_node(lst):
    numbers = []
    for i in range(len(lst)):
        name = lst[i].name.split('@')[0]
        numbers.append(int(name))
        max = numbers[0]
        for i in numbers:
            if max < i:
                max = i
    return(max)
    
def get_terminal_leaves(tree,post_order):
    df_t = []
    for clade in tree.get_nonterminals(order = post_order):
        df_t.append(clade.get_terminals(order = post_order))
        
    return(df_t)
    
def finding_duplicates(final_bit_vector):
    duplicateFrequencies = {}
    for i in set(final_bit_vector):
        duplicateFrequencies[i] = final_bit_vector.count(i)
        
    
    count = 0
    for values in (duplicateFrequencies):
        if duplicateFrequencies[values] == 2:
            count += 1
    return count

def compare_two(v1, v2):
    
   flag = True
   for i in range(len(v1)):
           if (v1[i] != v2[i]):
               flag = False
    
   return flag
    



if __name__=="__main__":

           
    if(len(sys.argv) < 3):
        sys.stderr.write('Its required two newick files to run the code. ')
        sys.exit(1)

    
    ###-------------------------TREE1--------------------------
    tree1 = Phylo.read(sys.argv[1], 'newick')
    tree2 = Phylo.read(sys.argv[2], 'newick')
               
    
    names = tabulate_names_one(tree1)
                    
    df_t1 = []
                    
    df_t1 = get_terminal_leaves(tree1,"postorder")

    lst = get_look_up_table(tree1)
    ###-------------------------------------------------------
                    
                    
    ####--------------------------TREE2-----------------------
    ###reading tree 2
    
    names = tabulate_names_two(tree2,lst)

                    
    df_t2 = []
    
    
                    
    for clade in tree2.get_nonterminals(order = "postorder"):
        min = min_node(clade.get_terminals(order = "postorder"))
        max = max_node(clade.get_terminals(order = "postorder"))
        size = len(clade.get_terminals(order = "postorder"))
        if max - min + 1 == size:
            df_t2.append(clade.get_terminals(order = "postorder"))
                    
                    
    
    df_t2_1 = get_terminal_leaves(tree2,"postorder")
    
    
    
    last = len(df_t2) -1 
    numLeaves = len(df_t2[last])
                    
    lstForT2 = []
                    
    for i in range(len(df_t2)):
                
        vector = create_bit_vector(df_t2[i],numLeaves)
        lstForT2.append(''.join(str(x) for x in vector))
                        
                
    lstForT1 = []
                    
    for i in range(len(df_t1)):
                
        vector = create_bit_vector(df_t1[i],numLeaves)
        lstForT1.append(''.join(str(x) for x in vector))
        
    bit_vector = lstForT2 + lstForT1
    bit_vector.sort()

    count = finding_duplicates(bit_vector)
    print("RF-distance : ", (len(df_t1) + len(df_t2_1) - 2*count))


    t = Tree()
    t.populate(10, random_branches=True)
    print(t)
    



    
    
    
            
            
            

    
    