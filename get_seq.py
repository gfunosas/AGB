
def FASTA_iterator(fasta_filename):
	"""A Generator Function that reads a Fastafile. In each iteration, the function must return a tuple with the following format:
	(identifier, sequence). Function name: FASTA_iterator(fasta_filename)"""
	
	fasta_file = open(fasta_filename)
	sequence = ""
	fasta_tuple = ()

	for line in fasta_file:
		if line[0]==">":
			if len(sequence) > 100 and sequence[-50] == 'G' and sequence[-49] == 'T':
				fasta_tuple=(name.split('|')[0], sequence[-100:])
				yield fasta_tuple
				sequence=""
			name=line.strip(">\n")
		else:
			sequence+=line.strip()
	if len(sequence) > 100 and sequence[-50] == 'G' and sequence[-49] == 'T':
		fasta_tuple=(name.split('|')[0],sequence[-100:])
	yield fasta_tuple

	fasta_file.close()


output = open('DNA_seq.txt', 'w')
	
for element in FASTA_iterator('mart_export.txt'):
	output.write('%s \n' %str(element))

