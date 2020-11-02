#!/usr/bin/env python3
# Author: Colin Naughton 
# This was originally written Spring 2016 for BME160: Research Programming in the Life Sciences with Dr. David Bernick.
'''
PSEUDOCODE:

Import sequenceAnalysis
Read FASTA file
Iterate through FASTA headers/sequences.
	Iterate through the different reading frames of each sequence.
		Iterate through each codon of sequence in a particular reading frame.
			If codon = start codon(top strand).
				Save codon position to list.
			If codon = stop codon(top strand).
				Save gene(frame, start codon position, end codon position).
			If codon position searched at end of sequence.
				If a start codon(top strand) has been found since last gene.
					Save gene(frame, start codon position, end of sequence)

			If codon = start codon(bottom strand).
				Save codon position to list.
			If codon = stop codon(top strand).
				If a start codon has been found since last gene.
					Save gene(frame, start codon position(previous), previous start codon position).
				Declare new stop codon(current position).
			If codon position searched at end of sequence.
				If a start codon(bottom strand) has been found since last gene.
					Save gene(frame, stop codon position(previous), sequence length).
				Reset start/stop codon lists(bottom strand).
	Sort genes in sequence by length, then start position.
	Print header/sequence gene info for genes>100 nucleotides.


'''
'''
Uses imported sequenceAnalysis module to iterate FASTA headers and their complimentary sequences
from FASTA file read from directory, finding open reading frames in the sequence.
input: fastA file
outputs: FASTA header along with all open reading frames above 100 nucleotides within header's complimentary sequence. 
'''

import sequenceAnalysis

class ORFstrand:

	def __init__(self,s):
		'''Capture the DNA strand and caste it in upper case. Initialize list to contain gene info as tuples.'''
		self.nucString = s.upper()
		self.TupleList = []
		'''Declare stop codons to be searched for.'''
		self.stopCodonList = ['TAA', 'TGA', 'TAG']
		self.revStopCodonList = ['TTA', 'TCA', 'CTA']

	def findORFs(self):

		'''Iterate through all 3 reading frames of a sequence, finding all ORFs in the sequence.'''
		for k in range(0, 3): 
			startCodonList = [1] #Initialize position of start codon at beginning of sequence being analyzed.
			revStopCodon = [1] #Initialize position of the opposite strand's stop codon at beginning of sequence being analyzed.
			revStartCodonList = [] #Initialize empty list of reverse start codons.

			'''Iterate through all codons in a reading frame, finding start and stop
			codons for the top and bottom strand.'''
			for i in range(k, len(self.nucString), 3): 
				codon = self.nucString[i:i + 3]

				'''Top-strand gene search'''
				if codon == 'ATG': #Find position of start codon on top strand, add to list.
					startCodonList.append(i + 1)

				if codon in self.stopCodonList: #Find position of stop codon on top strand, save gene info.
						
					if len(startCodonList) > 0:
							
						self.TupleList.append((k + 1, startCodonList[0], i + 3)) #Save gene info if a start codon is also present.
						startCodonList = [] #Clear list of top-strand start codons.

				'''Determine if top-strand search has reached end of sequence. Save gene info using remaining codon info'''
				if i >= len(self.nucString) - 3:
					if len(startCodonList) > 0:

						if len(startCodonList) > 1:
							if startCodonList[1] == 1: #Stop codon not found since first codon of reading frame.
								self.TupleList.append((k + 1, startCodonList[1], i + 3)) #

						elif startCodonList[0] != 1: #Stop codon found since first codon of reading frame.
							self.TupleList.append((k + 1, startCodonList[0], i + 3))

				'''Bottom-strand gene search'''
				if codon == 'CAT': #Find position of start codon on bottom strand, add to list.
					revStartCodonList.append(i + 3)

				'''If stop codon found that doesn't start at position 1, save gene info using previous bottom-strand stop codon.'''
				if codon in self.revStopCodonList:  
					if len(revStartCodonList) > 0:
						self.TupleList.append(( -(((len(self.nucString) - (revStopCodon[-1] - 1)) % 3) + 1), revStopCodon[-1], revStartCodonList[-1]))
						revStartCodonList = [] #Clear start codon list for bottom strand.
					revStopCodon.append(i + 1) #Declare new stop codon position.

				'''Determine if top-strand search has reached end of sequence. Save gene info using remaining codon info'''
				if i >= len(self.nucString) - 3:
					if len(revStopCodon) > 1 :
						self.TupleList.append((-(((len(self.nucString)-(revStopCodon[-1]-1))%3)+1), revStopCodon[-1], len(self.nucString)))
					revStopCodon = [1] #Reset stop codon list for bottom strand.
					revStartCodonList = [] #Clear list of bottom-strand start codons.

		return self.TupleList


########################################################################################################################################

'''Create object of FASTA headers and their complimentary sequences from FASTA file read from directory.'''
thing = sequenceAnalysis.FastAreader('tass2.fa')

'''Iterate through FASTA header(s) and sequence(s), finding gene(s) in the ORF.''' 
for head, seq in thing.readFasta() :

	testSequence = ORFstrand(seq)
	ORFs = testSequence.findORFs()

	'''Sort found genes by length, then start position'''
	newTupleList = [] #Initialize new list for sorted genes	
	for gene in sorted(ORFs, key=lambda x: ((x[2] - x[1]), x[1]), reverse=True): #Sort by gene length, then first position of gene.
		if (gene[2] - gene[1] + 1) >= 100:
			newTupleList.append(gene)

	'''Print header w/ relevant genes from the sequence.'''
	print(head) #Print header.
	for gene in newTupleList: #Iterate through genes, printing frame, start, stop, and length.
		print('{:+d} {:>5d}..{:>5d} {:>5d}'.format((gene[0]), (gene[1]), (gene[2]), (gene[2] - gene[1] + 1)))

    