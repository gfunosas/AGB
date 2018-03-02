import sys
from itertools import islice

#### Put the model you want to run:
# --------------------------------

#from model1 import get_donor
#from model2 import get_donor
#from model3_2 import get_donor
from model4_2 import get_donor

# ---------------------------------

module = sys.modules[get_donor.__module__]
module = str(module).split()[1].strip("'")


if module == 'model1':
	output = open('model1.txt', 'w')
elif module == 'model2':
	output = open('model2.txt', 'w')
elif module == 'model3_2':
	output = open('model3.txt', 'w')
elif module == 'model4_2':
	output = open('model4.txt', 'w')

num_test = 10
index = 0

with open('prova.txt', 'r') as infile:
	total_seq = len(infile.readlines())
	len_test = int(total_seq/num_test)
	Ni = 0
	Nf = len_test

with open('prova.txt', 'r') as infile:
	output.write(' MODEL\tTP\tFP  \tTPR  \tFDR\n')
	for i in range(0, num_test):
		FP = 0
		TP = 0
		
		lines_gen = islice(infile, Ni, Nf, 1)
		for line in lines_gen:
			line = line.replace("(", "").replace(")", "").strip()
			line = line.split(",")
			DNA_seq = line[1].strip().strip("'")
			
			if 'N' in DNA_seq:
				pass
			else:
				donor_pos = float(get_donor(DNA_seq))
				if donor_pos == 51:
					TP += 1
				elif donor_pos == 0:
					pass
				else:
					FP +=1
		try:			
			TPR = float(TP)/(float(TP)+float(FP))
			FDR = float(FP)/(float(TP)+float(FP))
			output.write('%s\t%s\t%s\t%.5f\t%.5f\n' %(module, TP, FP, TPR, FDR))
		except ZeroDivisionError:
			TPR = '-'
			FDR = '-'
			output.write('%s\t%s\t%s\t%s\t%s\n' %(module, TP, FP, TPR, FDR))
		index += TPR
	print(index/10)
