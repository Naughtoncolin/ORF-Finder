#!/usr/bin/env python3
# Authors: Colin Naughton & Samuel Pinsky-Dickson
# This code was written in Spring 2016 during BME160: Research Programming in the Life Sciences and adapted from code written by Dr. David Bernick.

"""
Modules:
1)NucParams: Models strand of DNA/RNA.

2)FastAreader: Reads file containing one or more FASTA
formatted sequences.

3)ProteinParam: Models an amino acid sequence.
"""

class NucParams:

    """
    Models a strand of DNA/RNA.

    Keyword arguements:
    polyNuc = String representing the nucleotides of a strand of DNA/RNA

    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
 
    def __init__ (self, polyNuc):
        """Capture the DNA/RNA strand and caste it in upper case."""
        l = ''.join(polyNuc).split()
        self.nucString = ''.join(l).upper()

        """Build list of codons comprising the DNA/RNA strand."""
        #self.nucCodons = [self.nucString[i:i+3] for i in range(0, len(self.nucString), 3)]
        self.nucCodons = {codon:0 for codon in NucParams.rnaCodonTable.keys()}
        self.aaCount = {aa:0 for aa in NucParams.rnaCodonTable.values()}
        self.nucCount = {'A':0,'T':0,'C':0,'G':0,'U':0, 'N':0}

    def addSequence (self, thisSequence):
        newSeq = thisSequence.replace('T', 'U')
        for i in range(0, len(thisSequence), 3):
            codon = newSeq[i:i+3]
            if codon in NucParams.rnaCodonTable.keys(): # change to in self.nucCodons
                self.nucCodons[codon] += 1
                self.aaCount[NucParams.rnaCodonTable[codon]] += 1
        for nuc in thisSequence:
            if nuc in self.nucCount.keys():
                self.nucCount[nuc] += 1



    def aaComposition(self):
        return self.aaCount

    def nucComposition(self):
        return self.nucCount

    def codonComposition(self):
        return self.nucCodons

    def nucCount1(self):
        """Returns count of the total number of nucleic acids in 
        a given DNA/RNA strand."""
        nucTotal = 0
        for key in self.nucCount.keys():
            nucTotal += self.nucCount[key]
        return nucTotal

 
 
class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname
 
    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''
        
        with open(self.fname) as fileH:
            # initialize return containers
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
 
# presumed object instantiation and example usage
# myReader = FastAreader ('testTiny.fa');
# for head, seq in myReader.readFasta() :
#     print (head,seq)


class ProteinParam :

    """
    Models an amino acid sequence.

    Keyword arguements:
    protein = String representing amino acids of a polypeptide.

    """

# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
 
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated
    def __init__ (self, protein):
        l = ''.join(protein).split()
        self.protString = ''.join(l).upper()
        self.protString = list(self.protString)
        
        newList = []    #Is there a better way to do this???
        for aa in self.protString:
            if aa  in ProteinParam.aa2mw.keys():
               newList.append(aa)
        self.protString = ''.join(newList)

    def aaCount (self):
        """Counts the amino acids of a given polypeptide."""
        return len(self.protString)
 
    def pI (self):
        """Returns the isoelectric point of a given polypeptide."""
        
        return self._charge_()

 
    def aaComposition (self) :
        """Determines percent composition of each amino acid in a given polypeptide."""
        aaCount = dict()
        for aa in ProteinParam.aa2mw.keys():
            aaCount[aa] = self.protString.count(aa)
        return aaCount

    def _charge_ (self):
        """Determine charge of a polypeptide at a certain pH."""
        
        isoelectricPoint = 0
        isoelectricCharge = 1000000000
        
        for pH in range(1401):
            netCharge = 0
            for aa in ProteinParam.aa2chargePos.keys():
                netCharge += (self.protString.count(aa) * 10 ** ProteinParam.aa2chargePos[aa]) / (10 ** ProteinParam.aa2chargePos[aa] + 10 ** (pH/100))
            for aa in ProteinParam.aa2chargeNeg.keys():
                netCharge -= (self.protString.count(aa) * 10 ** (pH/100)) / (10 ** ProteinParam.aa2chargeNeg[aa] + 10 ** (pH/100))
            netCharge +=  (10 ** ProteinParam.aaNterm) / (10 ** ProteinParam.aaNterm + 10  ** (pH/100))
            netCharge -= (10 ** (pH/100)) / (10 ** ProteinParam.aa2chargeNeg[aa] + 10 ** (pH/100))  
            if abs(netCharge) < abs(isoelectricCharge):
                isoelectricCharge = netCharge
                isoelectricPoint = (pH/100)
        return isoelectricPoint

    def molarExtinction (self):
        """Determine molar extinction coefficient of a given polypeptide."""
        extintCoef = 0
        for aa in ProteinParam.aa2abs280.keys():
            extintCoef += self.protString.count(aa) * ProteinParam.aa2abs280[aa]
        return extintCoef
 
    def massExtinction (self):
        """Determine the mass extinction coefficient of a given polypeptide."""
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
 
    def molecularWeight (self):
        """Determines the molecular weight of a given polypeptide."""
        protMass = 0
        for aa in ProteinParam.aa2mw.keys():
            protMass += (ProteinParam.aa2mw[aa] * self.protString.count(aa))
        protMass -= ProteinParam.mwH2O * (myParamMaker.aaCount() - 1)
        return protMass

#Test: VLSPADKTNVKAAW
"""Prints info on given polypeptide."""
# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
'''for inString in sys.stdin :
    myParamMaker = ProteinParam(inString)
    myAAnumber = myParamMaker.aaCount()
    print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
    print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
    print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
    print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
    print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
    print ("Amino acid composition:")
    myAAcomposition = myParamMaker.aaComposition()
    keys = list(myAAcomposition.keys())
    keys.sort()
    if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
    for key in keys :
        print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
    '''
    #myParamMaker = NucParams(inString)
    # presumed object instantiation and example usage
#myReader = FastAreader ('testGenome.fa.txt');
#for head, seq in myReader.readFasta() :
#    print (head)
