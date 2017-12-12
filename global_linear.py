import numpy as np
import sys


# Define parameters

# For longest common subsequence
#score_matrix = np.zeros((4, 4), int)
#np.fill_diagonal(score_matrix, 2)
#gap = -1

# For unit cost edit distance - change to MIN
#score_matrix = np.full((4, 4), 1, dtype=int)
#np.fill_diagonal(score_matrix, 0)
#gap = 1

#MANUAL INPUT
score_matrix = np.array([[10, 2, 5, 2], [2, 10, 2, 5], [5, 2, 10, 2], [2, 5, 2, 10]])
gap = 0

to_num = {"A": 0, "C": 1, "G": 2, "T": 3}
to_base = {0: "A", 1: "T", 2: "C", 3: "G"}

print "\nSubstitution matrix: \n", score_matrix, "\n"
print "Gap penalty: ", gap, "\n"

S1 = str(sys.argv[1])
S2 = str(sys.argv[2])
#S1 = "GGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCGGCATTTTAAGATAAGTGAGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACTTAAGCGTATGCCCAGATCTGACTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTCTGTTTACGTATAAACAGAATCGCCTGGGTTCGC"
#S2 = "GGGCTAAAGGTTAGGGTCTTTCACACTAAAGAGTGGTGCGTATCGTGGCTAATGTACCGCTTCTGGTATCGTGGCTTACGGCCAGACCTACAAGTACTAGACCTGAGAACTAATCTTGTCGAGCCTTCCATTGAGGGTAATGGGAGAGAACATCGAGTCAGAAGTTATTCTTGTTTACGTAGAATCGCCTGGGTCCGC"
n = len(S1)
m = len(S2)
print(n)

def score(i, j):
	return score_matrix[to_num[S1[i-1]], to_num[S2[j-1]]]


# Fill up table

def table(n, m):
	T = np.empty((n + 1, m + 1), int)

	for i in range(0, n+1):
		for j in range(0, m+1):

			if i == 0 and j > 0:
				T[i, j] = T[i, j-1] + gap

			if i > 0 and j > 0:
				T[i, j] = max(T[i-1, j-1] + score(i, j), T[i, j-1] + gap, T[i-1, j] + gap)

			if i > 0 and j == 0:
				T[i, j] = T[i-1, j] + gap

			if i == 0 and j == 0:
				T[i, j] = 0
	return T



# def backtrack(i, j):

# 	if i > 0 and j > 0 and T[i, j] == T[i-1, j-1] + score(i, j):
# 		#O1.append(S1[i-1]), O2.append(S2[j-1])
# 		#print i,j
# 		#backtrack(i-1, j-1)
# 	 	s1, s2  = backtrack(i-1, j-1)
# 	 	return s1 + S1[i-1], s2 + S2[j-1]

# 	if i >= 0 and j > 0 and T[i, j] == T[i, j-1] + gap: #GAP IN S1 - LEFT DIR
# 		#O1.append("-"), O2.append(S2[j-1])
# 		#print i,j
# 		#backtrack(i, j-1)
#  		s1, s2 = backtrack(i, j-1)
#  		return s1 + "-", s2 + S2[j-1]

# 	if i > 0 and j >= 0 and T[i, j] == T[i-1, j] + gap: # GAP IN S2 - UP DIR
# 		#O1.append(S1[i-1]), O2.append("-")
# 		#print i,j
# 		#backtrack(i-1, j)
#  		s1, s2 = backtrack(i-1, j)
# 		return s1 + S1[i-1], s2 + "-"

#  	if i == 0 and j == 0:
#  		#print i,j
#  		#O1.append("*"), O2.append("+")
#  		return "", ""



T = table(n, m)
print "Score matrix: \n", T

print "\nThe score of an optimal alignment is:", T[n][m]

# O1 = []
# O2 = []
# alignment = backtrack(n, m)

# seq = range(0, len(alignment[0]), 50)


# c = 0
# for s in seq:
# 	if s < 198:

# 		print " "
# 		print alignment[0][s:seq[c+1]]
# 		print alignment[1][s:seq[c+1]]
# 		c +=1
# 	else:
# 		print " "
# 		print alignment[0][s::]
# 		print alignment[1][s::]



# O1 = "".join(O1[::-1]).split("*")
# O2 = "".join(O2[::-1]).split("+")

# # for o1, o2 in zip(O1, O2):
	
# # 	print o1
# # 	print o2
# # 	print ""

# # num_align = len(O1) -1

# # #print "Number of global alignmets with max/min score: ", num_align



# def backtrack_n(i, j):

# 	if i > 0 and j > 0 and T[i, j] == T[i-1, j-1] + score(i, j):
# 		n_d_1 = S1[i-1]
# 		n_d_2 = S2[j-1]
# 		backtrack(i-1, j-1)

# 	if i >= 0 and j > 0 and T[i, j] == T[i, j-1] + gap: #GAP IN S1 - LEFT DIR
# 		n_l_1 = "-"
# 		n_l_2 = S2[j-1]
# 		backtrack(i, j-1)

# 	if i > 0 and j >= 0 and T[i, j] == T[i-1, j] + gap: # GAP IN S2 - UP DIR
# 		n_u_1 = S1[i-1]
# 		n_u_2 = "-"
# 		backtrack(i-1, j)

#  	if i == 0 and j == 0:
#  		#print i,j
#  		#O1.append("*"), O2.append("+")
#  		return "", ""



