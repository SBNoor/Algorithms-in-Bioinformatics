import sys
import numpy as np

## A. Functions

# A.1. Pairwise cost from substitution matrix
def sub(i, j, S1, S2):
	to_num = {"A": 0, "C": 1, "G": 2, "T": 3}
	to_base = {0: "A", 1: "T", 2: "C", 3: "G"}
	return score_matrix[to_num[S1[i-1]], to_num[S2[j-1]]]

# A.2. SP-cost of a column of 3 symbols
def SP(p1, p2, p3, i, j, k):
	if p1 == 1 and p2 == 1 and p3 == 1:
		return sub(i, j, seq1, seq2) + sub(j, k, seq2, seq3) + sub(i, k, seq1, seq3)
	if p1 == 1 and p2 == 1 and p3 == 0:
		return sub(i, j, seq1, seq2) + gap + gap
	if p1 == 1 and p2 == 0 and p3 == 1:
		return sub(i, k, seq1, seq3) + gap + gap
	if p1 == 0 and p2 == 1 and p3 == 1:
		return sub(j, k, seq2, seq3) + gap + gap
	else:
		return gap + gap

# A.3. Fill up dynamic programming table of n**3 dimensions
def dynamic(n1, n2, n3):

	T = np.empty((n1+1, n2+1, n3+1), int)

	for i in range(n1+1):
		for j in range(n2+1):
			for k in range(n3+1):
				v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("Inf")

				if i==0 and j==0 and k==0:
					v0 = 0
				if i>0 and j>0 and k>0:
					v1 = T[i-1][j-1][k-1] + SP(1, 1, 1, i, j, k)
				if i>0 and j>0 and k>=0:
					v2 = T[i-1][j-1][k] + SP(1, 1, 0, i, j, k)
				if i>0 and j>=0 and k>0:
					v3 = T[i-1][j][k-1] + SP(1, 0, 1, i, j, k)
				if i>=0 and j>0 and k>0:
					v4 = T[i][j-1][k-1] + SP(0, 1, 1, i, j, k)
				if i>0 and j>=0 and k>=0:
					v5 = T[i-1][j][k] + SP(1, 0, 0, i, j, k)
				if i>=0 and j>0 and k>=0:
					v6 = T[i][j-1][k] + SP(0, 1, 0, i, j, k)
				if i>=0 and j>=0 and k>0:
					v7 = T[i][j][k-1] + SP(0, 0, 1, i, j, k)

				T[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)
	return T

# A.4. Backtracking through the dynamic programming table
def backtrack(i, j, k, T):

	if i>0 and j>0 and k>0 and T[i][j][k] == T[i-1][j-1][k-1] + SP(1, 1, 1, i, j, k):
		s1, s2, s3 = backtrack(i-1, j-1, k-1, T)
		return s1 + seq1[i-1], s2 + seq2[j-1], s3 + seq3[k-1]

	if i>0 and j>0 and k>=0 and T[i][j][k] == T[i-1][j-1][k] + SP(1, 1, 0, i, j, k):
		s1, s2, s3 = backtrack(i-1, j-1, k, T)
		return s1 + seq1[i-1], s2 + seq2[j-1], s3 + "-"

	if i>0 and j>=0 and k>0 and T[i][j][k] == T[i-1][j][k-1] + SP(1, 0, 1, i, j, k):
		s1, s2, s3 = backtrack(i-1, j, k-1, T)
		return s1 + seq1[i-1], s2 + "-", s3 + seq3[k-1]

	if i>=0 and j>0 and k>0 and T[i][j][k] == T[i][j-1][k-1] + SP(0, 1, 1, i, j, k):
		s1, s2, s3 = backtrack(i, j-1, k-1, T)
		return s1 + "-", s2 + seq2[j-1], s3 + seq3[k-1]

	if i>0 and j>=0 and k>=0 and T[i][j][k] == T[i-1][j][k] + SP(1, 0, 0, i, j, k):
		s1, s2, s3 = backtrack(i-1, j, k, T)
		return s1 + seq1[i-1], s2 + "-", s3 + "-"

	if i>=0 and j>0 and k>=0 and T[i][j][k] == T[i][j-1][k] + SP(0, 1, 0, i, j, k):
		s1, s2, s3 = backtrack(i, j-1, k, T)
		return s1 + "-", s2 + seq2[j-1], s3 + "-"

	if i>=0 and j>=0 and k>0 and T[i][j][k] == T[i][j][k-1] + SP(0, 0, 1, i, j, k):
		s1, s2, s3 = backtrack(i, j, k-1, T)
		return s1 + "-", s2 + "-", s3 + seq3[k-1]

	if i==0 and j==0 and k==0:
		return "", "", ""

# A.5. Fasta file read
def reading_fasta_file(filename):
	f = open(filename,"r")
	line = f.readlines()[2:]
	f.close()
	sequences = []

	c = -1
	for i in range(len(line)):
		if line[i][0] == ">":
			name = line[i][1:].strip()
			sequences.append([name, ""])
			c += 1
		else:
		 	sequences[c][1] += line[i].rstrip()
	return sequences

# A.6. Fasta alignment write
def writing_fasta_file(in_seqs, out_seqs, filename):
	f = open(filename, "w")
	f.write("Optimal alignment for sequences in {0}\nScore of alignment = {1}\n\n".format(file1, score))
	n = 70
	for s in range(len(in_seqs)):
		name = in_seqs[s][0]
		out = out_seqs[s]
		if len(out) < n:
			f.write(">{0}\n{1}\n".format(name, out))
		else:
			out_s = [out[i:i+n] for i in range(0, len(out), n)]
			f.write(">{0}\n".format(name))
			for o in out_s:
				f.write("{0}\n".format(o))

## B. Program body

if __name__ == '__main__':
	score_matrix = np.array([[0, 5, 2, 5], [5, 0, 5, 2], [2, 5, 0, 5], [5, 2, 5, 0]])
	gap = 5

	# Set inputs
	file1 = sys.argv[1]
	in_seqs = reading_fasta_file(file1) #List of lists
	seq1, seq2, seq3 = [s[1] for s in in_seqs]
	n1 = len(seq1)
	n2 = len(seq2)
	n3 = len(seq3)

	# Built dynamic programming table
	T_dyn = dynamic(n1, n2, n3)
	score = T_dyn[n1][n2][n3]
	# Backtrack 1 optimal alignment
	out_seqs = backtrack(n1, n2, n3, T_dyn)

	# Write output
	writing_fasta_file(in_seqs, out_seqs, "{0}_MSA_exact.txt".format(file1[:-4]))
