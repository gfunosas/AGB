

hmm =	{'states':["begin","exon","donor","intron"],
		'emp':	[[0.00,0.00,0.00,0.00],	# begin is a silent state
				[0.25,0.25,0.25,0.25],	# exon nucleotides are equiprobable
				[0.05,0.00,0.95,0.00],	# donors sites are mostly G
				[0.40,0.10,0.10,0.40]],	# introns are A-Trich,
		'trp': 	[[0.0,1.0,0.0,0.0],	# out going from begin
				[0.0,0.9,0.1,0.0],	# out going from exon
				[0.0,0.0,0.0,1.0],	#out going from donor
				[0.0,0.0,0.0,1.0]],	#out going from intron
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
				possibilities = [1,2]
				for i in possibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]
			elif final_hmm[0] == 'donor':
				my_dic[Matrix[2][0]] = Matrix[2][col]
			elif final_hmm[0] == 'intron':
				possibilities = [3,4]
				for i in possibilities:
					my_dic[Matrix[i][0]] = Matrix[i][col]
			else: 
				print('not known state')		
				
		return (my_dic)
	
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
	print('\nThe donor position is %i \n' %donor_position)
else: 
	print('\nThere is no donor position \n')

