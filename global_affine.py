import numpy as np
import re
import time
import sys


###The following function returns the index of a base(in question) in the scoring matrix
def matchorMismatch(x):
	x = x.lower()
	i = 0
	if x == 'a':
		i = dictionary_bases['A']
		return i
	elif x == 'c':
		i = dictionary_bases['C']
	elif x == 'g':
		i = dictionary_bases['G']
	elif x == 't':
		i = dictionary_bases['T']

	return i

###--------------------------------- END OF FUNCTION----------------------------------------------------

###The following function fills a particular cell in the deletion matrix when called from matrix_fill function
def filling_del_matrix(i,j,dynamic_matrix,deletion_matrix):
	alpha = 5
	beta = 5
	score1 = score2 = None
	if i > 0 and j >= 0:
		score1 = dynamic_matrix[i-1][j] + (alpha + beta)
	if i > 1 and j >= 0:
		score2 = deletion_matrix[i-1][j] + alpha

	list = []
	list.append(score1)
	list.append(score2)

	list = [i for i in list if i is not None]

	deletion_matrix[i][j] = min(list)
	return(min(list))

###--------------------------------- END OF FUNCTION ----------------------------------------------------
	
###The following function fills the insertion matrix when called from the matrix_fill method
def filling_insertion_matrix(i,j,dynamic_matrix,insertion_matrix):
	alpha = 5
	beta = 5

	score1 = score2 = None

	if i >= 0 and j > 0:
		score1 = dynamic_matrix[i][j-1] + (alpha + beta)
	if i >= 0 and j > 1:
		score2 = insertion_matrix[i][j-1] + alpha

	list = []

	list.append(score1)
	list.append(score2)

	list = [i for i in list if i is not None]

	insertion_matrix[i][j] = min(list)
	return(min(list))

###--------------------------------- END OF FUNCTION -----------------------------------------------------

### The following function fills the dynamic programming table
def matrix_fill(dynamic_matrix,deletion_matrix,insertion_matrix,sequence_one,sequence_two,len_seq_one,len_seq_two,dictionary_bases,scoring_matrix):
	dynamic_matrix[0][0] = 0
	insertion_matrix[0][0] = 0
	deletion_matrix[0][0] = 0

	score_diag = score_left = score_up = None

	for i in range(0,len_seq_one):
		for j in range(0, len_seq_two):
			if i > 0 and j > 0:

				row = matchorMismatch(sequence_one[i-1])
				col = matchorMismatch(sequence_two[j-1])

				score_diag = dynamic_matrix[i-1][j-1] + scoring_matrix[row][col]

			if i > 0 and j >= 0:
				score_up = filling_del_matrix(i,j,dynamic_matrix,deletion_matrix)
			if i >= 0 and j > 0:
				score_left = filling_insertion_matrix(i,j,dynamic_matrix,insertion_matrix)

			list = []
			list.append(score_diag)
			list.append(score_left)
			list.append(score_up)

			list = [i for i in list if i is not None]

			if not list:
				dynamic_matrix[i][j] = 0

			else:
				dynamic_matrix[i][j] = min(list)

	return(dynamic_matrix)

###------------------------------- END OF FUNCTION ---------------------------------------------

def trace_back(sequence_one_list,sequence_two_list,dynamic_matrix):
	i = len(sequence_one_list) 
	j = len(sequence_two_list) 

	seqOne = []
	seqTwo = []

	alpha = 5
	beta = 5

	while(i>0 or j>0):
		row = matchorMismatch(sequence_one[i-1])
		col = matchorMismatch(sequence_two[j-1])

		if (i > 0 and j > 0) and (dynamic_matrix[i][j] == dynamic_matrix[i-1][j-1] + scoring_matrix[row][col]):
			seqOne += sequence_one[i-1]

			seqTwo += sequence_two[j-1]

			i = i-1
			j = j-1
		else:
			k = 1
			while True:
				if (i>=k) and (dynamic_matrix[i][j] == dynamic_matrix[i-k][j] + (alpha + beta*k)):
					substring = sequence_one_list[(i-k):i]
					substring1 = substring[::-1]
					seqOne.extend(substring1)
					x = '-'*(k)
					seqTwo += x
					i = i - k
					break

				elif(j >= k) and (dynamic_matrix[i][j] == dynamic_matrix[i][j-k] + (alpha + beta*k)):
					substring = sequence_two_list[(j-k):j]
					substring1 = substring[::-1]
					seqTwo.extend(substring1)
					x = '-'*(k)
					seqOne += x
					j = j - k 
					break

				else:
					k = k + 1


	print(">Seq1")
	str1 = ''.join(seqOne)
	print(str1[::-1])
	print("\n>Seq2")
	str2 = ''.join(seqTwo)
	print(str2[::-1])

##### ------------------------------- END OF FUNCTION ------------------------------------------

def reading_fasta_file(filename):
	f = open(filename,'r')
	line = f.readlines()[1:]

	f.close()

	sequence = ''
	for i in range(len(line)):
		sequence = sequence + line[i].rstrip()
	return sequence

### ---------------------------------END OF FUNCTION ------------------------------------------

def extracting_scoring_matrix(filename,dictionary_bases):
	size_of_scoring_matrix = 0
	counter = 0
	row_matrix = 0

	with open(filename) as f:
		line = f.readline()
		size_of_scoring_matrix = int(line.strip())
		scoring_matrix = np.zeros([size_of_scoring_matrix,size_of_scoring_matrix])
		line = f.readline()

		while line:
			line = line.strip()
			x = re.split('\s+',line)
			key = x[0]
			dictionary_bases[key] = counter
			counter += 1
			for i in range(1,len(x)):
				scoring_matrix[row_matrix][i-1] = int(x[i])

			row_matrix += 1

			line = f.readline()

	return scoring_matrix,dictionary_bases

### --------------------------------- END OF FUNCTION ---------------------------------------

if __name__ == '__main__':

	###Taking user input
	file_one = sys.argv[1]
	file_two = sys.argv[2]
	file_three = sys.argv[3]

	dictionary_bases = {}
	scoring_matrix = extracting_scoring_matrix(file_three,dictionary_bases)[0]
	dictionary_bases = extracting_scoring_matrix(file_three,dictionary_bases)[1]


	sequence_one = reading_fasta_file(file_one)
	sequence_two = reading_fasta_file(file_two)

	###converting strings to list
	sequence_one_list = list(sequence_one)
	sequence_two_list = list(sequence_two)

	len_seq_one = len(sequence_one) + 1
	len_seq_two = len(sequence_two) + 1

	###initializing the 3 matrices
	dynamic_matrix = np.zeros([len_seq_one,len_seq_two])
	deletion_matrix = np.zeros([len_seq_one,len_seq_two])
	insertion_matrix = np.zeros([len_seq_one,len_seq_two])

	#start_time = time.clock()
	dynamic_matrix = matrix_fill(dynamic_matrix,deletion_matrix,insertion_matrix,sequence_one,sequence_two,len_seq_one,len_seq_two,dictionary_bases,scoring_matrix)
	print("Optimal Cost : ", dynamic_matrix[len_seq_one-1][len_seq_two-1])
	#time_elapsed = time.clock() - start_time
	#print(time_elapsed)

	#start_time = time.clock()
	trace_back(sequence_one_list,sequence_two_list,dynamic_matrix)
	#print("traceback : ", time.clock()-start_time)


