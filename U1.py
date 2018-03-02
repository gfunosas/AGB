import math
import pprint
from itertools import islice

U1 = 'CAGGTAAGT'
nucleotides = 'ACGT'
testing = {}
num_test = 10

with open('DNA_seq.txt', 'r') as infile:
	total_seq = len(infile.readlines())
	len_test = int(total_seq/num_test)
	Ni = 0
	Nf = len_test

for x in range(0, num_test-1):
	with open('DNA_seq.txt', 'r') as infile:
		FP, TP, num_seq  = 0,0,0
		dic, S, I = {}, {}, {}
		model = ('model%s' %x)
		testing[model] = islice(infile, Ni, Nf, 1)
		for line in infile:
			if line not in testing[model]:
				num_seq += 1
				sequence = line.split(",")[1].replace("'", "").replace(")", "").strip(" ")
				positions = sequence[42:60]
				
				for i in range(len(positions)):
					if i not in dic:
						dic[i] = {}
					else:
						for nucleotide in nucleotides:
							if positions[i] == nucleotide:
								if nucleotide not in dic[i]:
									dic[i][nucleotide] = 1
								else:
									dic[i][nucleotide] += 1

		for dic_i in dic:
			for key in dic[dic_i]:
				dic[dic_i][key] = round(dic[dic_i][key]/(num_seq-1), 4)
								
		for j in range (len(positions)):
			S[j] = 0
			for nucleotide in nucleotides:
				if nucleotide not in dic[j]:
					S[j] += 0
				elif dic[i][nucleotide] == 1:
					S[j] = 0
					break
				else:
					S[j] -= math.log(dic[j][nucleotide],2)*dic[j][nucleotide]
		for key in S:
			I[key] = round(math.log(4, 2) - S[key], 3)
		
		Ni += len_test
		Nf += len_test
		
		print('\n')
		pprint.pprint(dic)
		pprint.pprint(I)
