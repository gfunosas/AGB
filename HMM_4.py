
hmm =	{'states':["begin","exon", "S-1", "S-2", "S-3", "S-4", "S-5","ex-48", "ex-49", "ex-50", "donor", "int-T", "int-53", "int-54", "int-55", "int-56", "intron", "T-st"],
	'emp':	[[0.00,0.00,0.00,0.00],	# begin is a silent state
			[0.25,0.25,0.25,0.25],	# exon nucleotides are equiprobable
			[0.40,0.10,0.40,0.10], # S-1
			[0.40,0.10,0.40,0.10], # S-2
			[0.38,0.12,0.38,0.12], # S-3
			[0.38,0.12,0.38,0.12], # S-4
			[0.38,0.12,0.38,0.12], # S-5
			[0.32,0.36,0.19,0.13],  # ex-48
			[0.62,0.11,0.13,0.14],	# ex-49
			[0.10,0.03,0.80,0.07], 	# ex-50
			[0.00,0.00,1.00,0.00],	# donor-G
			[0.00,0.00,0.00,1.00],  # int-T
			[0.55,0.04,0.37,0.04],	# int-53
			[0.67,0.08,0.13,0.12],	# int-54
			[0.09,0.06,0.77,0.08],	# int-55
			[0.18,0.16,0.20,0.46],	# int-56
			[0.40,0.10,0.10,0.40],	# intron
			[0.06,0.06,0.06,0.82]],	# T-st
	'trp': 	[[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],	# out going from begin
			[0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],  # out going from exon
			[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],	# out going from S-1
			[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from S-2
			[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from S-3
			[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from S-4
			[0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from S-5
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from ex-48
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],  # out going from ex-49
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from ex-50
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0], # out going from donor-G
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0], # out going from int-T
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0], # out going from int-53
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0], # out going from int-54
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0], # out going from int-55
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0], # out going from int-56
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.8,0.2], # out going from intron
			[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.2,0.8]], # out going from T-state
	'initstate': 0
	}
		
import math, sys

class IncorrectSequenceLetter(ValueError):
		def __init__ (self, letter, class_name):
			self.letter = letter
			self.class_name = class_name

		def __str__ (self):
			return "The sequence item %s is not found in the alphabet of class %s" %(self.letter, self.class_name)
			
class DNASequence(object):

	alphabet ='GATC'
	
	def __init__(self, sequence):

		self.__sequence = sequence
		
		for letter in self.__sequence:
			if letter not in self.alphabet:
				e = IncorrectSequenceLetter(letter, self.__class__.__name__)
				raise e
				
	def get_sequence(self):
		return self.__sequence

DNA_seq = DNASequence(sys.argv[1])
seq_input = '-' + DNA_seq.get_sequence()

Matrix = [[0 for x in range(0, len(seq_input))] for y in range(0, len(hmm['states'])+1)]

Matrix[0] = list(seq_input)

col1 = hmm['states']
col1.insert(0, 'query')

for row in range (0, len(Matrix)): 
	Matrix[row].insert(0, col1[row])

def calculate_values (row,col):
	
	first_col = False
	
	if Matrix[0][col] == 'A':
		nuc = 0
	elif Matrix[0][col] == 'C':
		nuc = 1
	elif Matrix[0][col] == 'G':
		nuc = 2
	elif Matrix[0][col] == 'T':
		nuc = 3
	else: 
		first_col = True
	
	def f_max (k):
		
		if Matrix[k][col-1] == '-inf':
			value = None
		else: 
			try:
				value = float(Matrix[k][col-1]) + math.log(hmm['trp'][k-1][row])
			except:
				value = None
		return value
	
	if first_col == True:
		if row == 1: 
			result = hmm['initstate']
		else:
			result = '-inf'
			
	else:
		row = row-1
		options = [f_max(k) for k in range(1, len(hmm['states'])) if f_max(k) != None]
		
		try:
			result = math.log(hmm['emp'][row][nuc]) + float(max(options))
		except ValueError:
			result = '-inf'
		
	return round(float(result),2)

for i in range(1, len(Matrix[row])):
	for row in range(1, len(Matrix)):
		Matrix[row][i] = calculate_values(row,i)

def get_hmm ():
	final_hmm = []	
	i = 1
		
	def get_dic (col):
		my_dic = {}		
		
		if col == -1:
			for i in range(1, len(Matrix)):
				my_dic[Matrix[i][0]] = Matrix[i][col]
		else:
			if final_hmm[0] == 'exon':
				posibilities = [1,7]
				for i in posibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]					
			elif final_hmm[0] == 'S-1':
				my_dic[Matrix[2][0]] = Matrix[2][col]
			elif final_hmm[0] == 'S-2':
				my_dic[Matrix[3][0]] = Matrix[3][col]
			elif final_hmm[0] == 'S-3':
				my_dic[Matrix[4][0]] = Matrix[4][col]
			elif final_hmm[0] == 'S-4':
				posibilities = [5,7]
				for i in posibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]
			elif final_hmm[0] == 'S-5':
				my_dic[Matrix[6][0]] = Matrix[6][col]
			elif final_hmm[0] == 'ex-48':
				my_dic[Matrix[2][0]] = Matrix[2][col]
			elif final_hmm[0] == 'ex-49':
				my_dic[Matrix[8][0]] = Matrix[8][col]
			elif final_hmm[0] == 'ex-50':
				my_dic[Matrix[9][0]] = Matrix[9][col]
			elif final_hmm[0] == 'donor':
				my_dic[Matrix[10][0]] = Matrix[10][col]
			elif final_hmm[0] == 'int-T':
				my_dic[Matrix[11][0]] = Matrix[11][col]
			elif final_hmm[0] == 'int-53':
				my_dic[Matrix[12][0]] = Matrix[12][col]
			elif final_hmm[0] == 'int-54':
				my_dic[Matrix[13][0]] = Matrix[13][col]
			elif final_hmm[0] == 'int-55':
				my_dic[Matrix[14][0]] = Matrix[14][col]
			elif final_hmm[0] == 'int-56':
				my_dic[Matrix[15][0]] = Matrix[15][col]
			elif final_hmm[0] == 'intron':
				possibilities = [16,17,18]
				for i in possibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]
			elif final_hmm[0] == 'T-st':
				possibilities = [17,18]
				for i in possibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]
			else:
				print('Not known state')

				
		return my_dic
	
	while i <= len(seq_input)-1:
		col = -i
		my_dic = get_dic(col)
		min_threshold = float('-inf')
						
		for key, value in my_dic.items():
			if str(value) == '-inf':
				pass
			if value > min_threshold:
				min_threshold = value
				state = key
			
		final_hmm.insert(0, state)
		i += 1
		
	try:
		donor = final_hmm.index('donor') +1
	except:
		donor = 0
		
	final_hmm.insert(0, '-')
	final_hmm.insert(0, 'HMM')
	return (final_hmm, donor)
		

for row in Matrix:
	print('%s \t %s' %(row[0], '  '.join('%07s' % row[i] for i in range(1,len(row)))))

print('\n%s' %('  '.join('%07s' % state for state in get_hmm()[0])))

donor_position = get_hmm()[1]
if donor_position != 0:
	print('The donor position is %i' %donor_position)
else: 
	print('There is no donor position')

