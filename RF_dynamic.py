import sys
from ete3 import Tree
import time


def create_dic(list_order):
	dic = {}
	for a, b in enumerate(list_order):
		dic[b] = a
	return dic

def node_size(node):
	out = []
	for n in node.iter_search_nodes():
		if n.is_leaf():
			out.append(n)
	return len(out)

def count_int(node):
	out = []
	for n in node.iter_search_nodes():
		if not n.is_leaf():
			out.append(n)
	return len(out)

def initialize(tree1, tree2):
	nodes1 = [node for node in tree1.iter_search_nodes()]
	indexes1 = create_dic(nodes1)
	sizes1 = [node_size(n) for n in nodes1]

	nodes2 = [node for node in tree2.iter_search_nodes()]
	indexes2 = create_dic(nodes2)
	sizes2 = [node_size(n) for n in nodes2]

	T = []
	for e in range(len(nodes1)):
		T.append([None for k in range(len(nodes2))])

	return T, nodes1, indexes1, sizes1, nodes2, indexes2, sizes2

def size(e1, e2):

	if T[e1][e2] == None:
		# Case 1: both are end leaf
		if n1[e1].is_leaf() and n2[e2].is_leaf():
			if n1[e1].name == n2[e2].name:
				value = 1
			else:
				value = 0
		# Case 2: e1 is end leaf, but e2 is internal
		elif n1[e1].is_leaf():
			value = 0
			for child in n2[e2].children:
				value += size(e1, i2[child])
		# Case 3: e2 is end leaf, but e1 is internal
		elif n2[e2].is_leaf():
			value = 0
			for child in n1[e1].children:
				value += size(i1[child], e2)
		# Case 4: both are internal nodes
		else:
			value = 0
			for child1 in n1[e1].children:
				for child2 in n2[e2].children:
					value += size(i1[child1], i2[child2])
		T[e1][e2] = value

	return T[e1][e2]

def RFd(tree1, tree2, T):
	shared = 0
	for node1 in tree1.iter_search_nodes():
		for node2 in tree2.iter_search_nodes():
			if not node1.is_leaf() and not node2.is_leaf():
				if T[i1[node1]][i2[node2]] != None:
					intersection = T[i1[node1]][i2[node2]]
					if intersection == s1[i1[node1]] and intersection == s2[i2[node2]]:
						shared += 1
				else:
					intersection = size(i1[node1], i2[node2])
					if intersection == s1[i1[node1]] and intersection == s2[i2[node2]]:
						shared += 1

	return count_int(n1[0]) + count_int(n2[0]) - 2*shared

if __name__ == '__main__':

	#tree1 = Tree("tree1.new")
	#tree2 = Tree("tree2.new")

	n = int(sys.argv[1])
	tree1 = Tree()
	tree1.populate(n, random_branches=True)
	tree2 = Tree()
	tree2.populate(n, random_branches=True)

	T, n1, i1, s1, n2, i2, s2 = initialize(tree1, tree2)

	#Fill up dynamic programming table worst-case O(n2)
	for i in range(len(n1)):
		for j in range(len(n2)):
			size(i,j)

	#Calculate RF distance - n_non_trival_splits squared
	distance = RFd(tree1, tree2, T)
	
	print("RF distance: {0}".format(distance))








