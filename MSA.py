import numpy
from itertools import combinations
import sys

def matchorMismatch(x):
    x = x.lower()
    i = 0
    if x == 'a':
        i = dictionary_bases['A']
    elif x == 'c':
        i = dictionary_bases['C']
    elif x == 'g':
        i = dictionary_bases['G']
    elif x == 't':
        i = dictionary_bases['T']
    return i
 
    
    
            
def sp_score(M, score_matrix):
    
    sumScore = 0
    
    for  elem in M:
        elem = list(elem)
        x = list(combinations(elem, 2))
        
        for item in x:
            if item[0] == '-' and item[1]== '-':
                sumScore += 0
            elif (item[0] != '-' and item[1] == '-') or (item[0] == '-' and item[1] != '-'):
                sumScore+= 5
            else:
                a = matchorMismatch(item[0])
                b= matchorMismatch(item[1])
                sumScore+= score_matrix[a][b]
                
        
    return sumScore
    

def matrix_fill(sequence_one, sequence_two):
	gap_cost = 5


	score_diag = score_left = score_up = 0

	len_seq_one = len(sequence_one) + 1
	len_seq_two = len(sequence_two) + 1

	array = numpy.zeros([len_seq_one,len_seq_two])

	array[0][0] = 0

	for i in range(1, len_seq_one):
		array[i][0] = gap_cost * i
		#print(array[i][0])

	for j in range(1, len_seq_two):
		array[0][j] = gap_cost * j

	for i in range(1,len_seq_one):
		for j in range(1, len_seq_two):
			row = matchorMismatch(sequence_one[i-1])
			column = matchorMismatch(sequence_two[j-1])
			score_diag = array[i-1][j-1] + scoring_matrix[row][column]
			score_left = array[i][j-1] + gap_cost
			score_up = array[i-1][j] + gap_cost

			array[i][j] = min(score_diag,score_left,score_up)

	return (traceback(array,sequence_one,sequence_two,len_seq_one,len_seq_two))

def traceback(array,sequence_one,sequence_two,len_seq_one,len_seq_two):
	gap_cost = 5

	

	seqOne = ""
	seqTwo = ""

	i = len_seq_one - 1
	j = len_seq_two - 1

	while i > 0 and j > 0:
		score_current = array[i][j]
		score_diagonal = array[i-1][j-1]
		score_up = array[i][j-1]
		score_left = array[i-1][j]

		row = matchorMismatch(sequence_one[i-1])
		column = matchorMismatch(sequence_two[j-1])

		if score_current == score_diagonal + scoring_matrix[row][column]:
			seqOne += sequence_one[i-1]
			seqTwo += sequence_two[j-1]
			i -= 1
			j -= 1

		elif score_current == score_up + gap_cost:
			seqOne += "-"
			seqTwo += sequence_two[j-1]
			j -= 1

		elif score_current == score_left + gap_cost:
			seqOne += sequence_one[i-1]
			seqTwo += "-"
			i -= 1

	while i > 0:
		seqOne += sequence_one[i-1]
		seqTwo += "-"
		i -= 1

	while j > 0:
		seqOne += "-"
		seqTwo += sequence_two[j-1]
		j -= 1


	return((seqOne[::-1]),(seqTwo[::-1]))

def findingCentralSequence(matrix_to_find_Sc,lst_of_sequences):
	### summing over all columns in each row

	sum = max = index = 0
	for i in range(len(matrix_to_find_Sc)):
		for j in range(len(matrix_to_find_Sc[i]) - 1):
			sum += matrix_to_find_Sc[i][j]

		matrix_to_find_Sc[i][j+1] = sum
		sum = 0

	col = len(lst_of_sequences) - 1

	for i in range(len(lst_of_sequences)):
		if max > matrix_to_find_Sc[i][col]:
			max = matrix_to_find_Sc[i][col]
			index = i

	return index

def minimumEditDistance(x,y,scoring_matrix):
	gap_penalty = 5
	D = numpy.zeros([len(x)+ 1,len(y) + 1],dtype = int)
	D[0,1:] = range(1,len(y) + 1)
	D[1:,0] = range(1,len(x) + 1)
	for i in range(1,(len(x) + 1)):
		for j in range(1,(len(y) + 1)):
			row = matchorMismatch(x[i-1])
			col = matchorMismatch(y[j-1])
			score_diag = D[i-1][j-1] + scoring_matrix[row][col]
			score_up = D[i-1][j] + gap_penalty
			score_left = D[i][j-1] + gap_penalty
			D[i][j] = min(score_diag,score_up,score_left)

	return D[len(x)][len(y)]

def pairwiseMatrix(lst_of_sequences,scoring_matrix):
	for i in range(len(lst_of_sequences)):
		for j in range(len(lst_of_sequences)):
			if i == j:
				matrix_to_find_Sc[i][j] = 0
			else:
				x = ''.join(lst_of_sequences[i])
				y = ''.join(lst_of_sequences[j])
				matrix_to_find_Sc[i][j] = minimumEditDistance(x,y,scoring_matrix)

	return(matrix_to_find_Sc)

def making_list(M):
    lst = []
    seq = ''.join(str(r) for r in M[0])
    for i in range(len(seq)):
        if(seq[i] == '-'):
            continue
        else:
            lst.append(i)
    return lst

def invariant_test(lst_alignments,M):
    lst_check = []

    counter = 0
    for i in range(1,len(M)):
        seq_one = M[0][1:len(M[0])-1]
        seq_two = M[i][1:len(M[i])-1]

        seqOne = ''
        seqTwo = ''
        for j in range(len(seq_one)):
            if seq_one[j] == '-' and seq_two[j] == '-':
                continue
            else:
                seqOne += seq_one[j]
                seqTwo += seq_two[j]

        lst_check.append((seqOne,seqTwo))
    
    
    for i in range(len(lst_alignments)):
        print("index : ",i)
        print(lst_alignments[i][0])
        print(lst_check[i][0])
        print()
        print(lst_alignments[i][1])
        print(lst_check[i][1])
        print()
            
        
    

    
def jump_function(M, Sc, Sj, gap):
    
    i = 0
    walk_M0 = 0
    while(i != gap):
        if(Sc[i] != '-'):
            while True:
                if(M[0][walk_M0] != '-'):
                    walk_M0+= 1
                    break
                else:
                    walk_M0 += 1

        i+=1
        
           
        
    for index in range(len(M[:len(M)-1])):

        M[index] = M[index][:walk_M0] + '-' + M[index][walk_M0:]


    M[len(M)-1] = M[len(M)-1][:walk_M0] + Sj[gap] + M[len(M)-1][walk_M0:]

        
    
def mergingAlignments(M,S1,Si):
    if len(M) == 0:
        Sc = matrix_fill(S1,Si)[0]
        Sj = matrix_fill(S1,Si)[1]
        
        lst_alignments.append((Sc,Sj))

        Sc = "#" + Sc + "#"
        Sj = "#" + Sj + "#"
        M.append(Sc)
        M.append(Sj)
        
    else:
        lst = making_list(M)
        Sc = matrix_fill(S1,Si)[0]
        Sj = matrix_fill(S1,Si)[1]
        
        lst_alignments.append((Sc,Sj))

        if '-' not in Sc:
            Sj = "#" + Sj + "#"
            Sk = ''
            for i in range(len(lst)-1):
                if((lst[i+1] - lst[i])-1 == 0):
                    Sk = Sk + Sj[i]
                elif((lst[i+1] - lst[i]) - 1 != 0):
                    Sk = Sk + Sj[i] + ("-"*((lst[i+1] - lst[i])-1))
            ##adding the last elemnet in Sk
            Sk += Sj[len(Sj) - 1]
            M.append(Sk)
            
        elif '-' in Sc:
            index = []
            Sc_new = ''
            Sj_new = ''
            ###Skipping the gaps in Sc and making a perfect alignment
            for i in range(len(Sc)):
                if Sc[i] == '-':
                    continue
                else:
                    Sj_new += Sj[i]
                    
                    
            Sj_new = '#' + Sj_new + "#"
            #print("Sj : ", Sj_new)
            Sk = ''
            for i in range(len(lst)-1):
                if((lst[i+1] - lst[i])-1 == 0):
                    Sk = Sk + Sj_new[i]
                elif((lst[i+1] - lst[i]) - 1 != 0):
                    Sk = Sk + Sj_new[i] + ("-"*((lst[i+1] - lst[i])-1))
            M.append(Sk+"#")
            
            
            Sc = '#' + Sc + '#'
            Sj = '#' + Sj + '#'
            
            for i  in range(len(Sc)):
                if(Sc[i] == '-'):
                    jump_function(M, Sc,Sj,i)
                

def reading_fasta_file(filename):
    try:
        inFile=open(filename,'r')
        #lines = f.readlines()
        #f.close()
    except IOError:
        print("This file does not exist")
        return
    
    headerList = []
    seqList = []
    currentSeq = ''
    for line in inFile:
        if line[0] == ">":
            headerList.append(line[1:].strip())
            if currentSeq != '':
                seqList.append(currentSeq)
            currentSeq = ''
        else:
            currentSeq += line.strip()
    seqList.append(currentSeq)
    
    return(seqList)
        



if __name__ == '__main__':
    dictionary_bases = {'A':0, 'C':1, 'G': 2, 'T': 3}
    
    filename = sys.argv[1]
    
    lst_of_sequences = reading_fasta_file(filename)
    
    scoring_matrix = [[0,5,2,5],
                      [5,0,5,2],
                      [2,5,0,5],
                      [5,2,5,0]]
    
    matrix_to_find_Sc = numpy.zeros([len(lst_of_sequences),(len(lst_of_sequences)+1)])
    matrix_to_find_Sc = pairwiseMatrix(lst_of_sequences,scoring_matrix)
    i = findingCentralSequence(matrix_to_find_Sc,lst_of_sequences)
    
    
    lst_alignments = []
    
    
    ###Adding the central sequence to the Matrix
    M = []
    for j in range(len(lst_of_sequences)):
        if j == i:
            continue
        else:
            mergingAlignments(M,lst_of_sequences[i],lst_of_sequences[j])
            
    ###writing fasta file
            
    file = open("output.fasta","w")
    
    for i in range(len(M)):
        file.write(">sequence ")
        file.write(str(i+1))
        file.write("\n")
        file.write(M[i][1:len(M[0])-1])
        file.write("\n")
        
    file.close()
    
    for i in range(len(M)):
        M[i] = list(M[i])
        
    print("Score : ",  sp_score(numpy.asarray(M).T.tolist(), scoring_matrix))


    








