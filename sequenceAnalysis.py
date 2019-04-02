#!/usr/bin/env python3
#Name: Patrick Stumps(pstumps)
#Group Members: none

import sys
import os
import operator
import math
#from findORFs import CommandLine
#from FastAreader import FastAreader

'''
sequenceAnalysis.py
This is my toolbox of bioninformatics methods. It contains the classes orfFinder, nucParams,fastAreader, and proteinparams. Proteinparams functions to 
identify specific parameters of a given string representing a peptide. FastAreader is a helper class to read headers and sequences from fasta files and yield
the sequence. NucParams is a genomic analysis module which identifies characteristics of a given gene. OrfFinder is functions to identify the longest open reading
frames of in given genomic sequence. For more information on each of these and methods please see their respective classes.
'''

class orfFinder:
    '''
    class orfFinder
    This is my interpretation of orfFinder. My workflow for designing this program started with attempting to find all start codons and stop codons in a given sequence
    and appending their index positions in the string to a respective start and stop codon list. However, many problems quickly arose when I found difficulty
    attempting to sort out the start codons found only after stop codons. However, my problems were solved by simply "saving" the start codon found in the first
    index position in a list, and then clearing that list once a stop codon was found using boolean flags.  This allowed precise excision of the ORFs, but also compromised the programs
    ability to find smaller nested ORFs. I'm sure theres some way to give the program that kind of functionality but I ran out of time. It does, however, support
    the functionality of multiple start and stop codons and also allows for only showing orfs of the given length specified by command line parameters.
    '''
    def __init__(self, longestGene, minGene, startCodons, stopCodons):
        '''
        OrfFinder init. Takes lG, minGene, startCodons, stopCodons from commandline in findOrfs.py and initialiazes everything.
        Important note: this does not accomodate for the longest Gene argument. However, it should accomodate for other start and stop codons, as well as the minimum gene given.
        '''

        self.longestGene = longestGene
        self.mGene = minGene
        self.startsC = startCodons
        self.stopsC = stopCodons
        self.orfDict = []
        #self.startsC = ["ATG", "TTG", "GTG"]
        #self.stopsC = ["TAA", "TGA", "TAG"]


    def findOrf(self, head, sequence):
        '''
        def findOrf. Takes the header and the sequence from fastAreader and finds orfs.
        My Pseudocode is as follows:
            for frames 1-3
                for i in range frame, length of sequence, step 3:
                    codon = slice string by 3
                    if codon = start
                        save index by appending to startlist
                        start found = true
                    if start codon found and stop codon found:
                        save orf from saved start position to index currently at
                        reset start list
                        found start = false
                    if stop codon found but no start
                        save orf from beginning

        The rest is the same for the reverse except for the reverse complement.
        This program allows for searching both the forward and reverse sequence at the same time. 
        '''
        #print("Finding ORFs")

        starts = [] #Start position list
        startFound = False #Found Start 
        codonFound = False #Found Stop

        rStarts = []  #Reverse start position list
        rStartFound = False #found Start
        rCodonFound = False #found Stop

        reverseTrans = str.maketrans("ATGC", "TACG") #translate
        revSequence = sequence[::-1].translate(reverseTrans) #backwards sequences

        for f in range(0,3): #Determines frame
            starts = [] #Saved start position
            startFound = False #found start
            codonFound = False #found stop

            rStarts = [] #saved reverse start
            rStartFound = False
            rCodonFound = False
            for i in range(f, len(sequence), 3): # search from frame to length of sequence
                codon = sequence[i:i+3] #codon slice
                rCodon = revSequence[i:i+3] #reverse codon slice
                #code for forward sequence begins
                if codon in self.startsC: #if start codon found, append
                    startFound = True 
                    starts.append(i)
                if codon in self.stopsC and startFound == True: # if stop codon found and start codon found
                    #print(head)
                    #print(starts)
                    if ((i+3) - (starts[0]+1))+1 >= self.mGene: #check of appropriate length
                        self.orfSaver((f%3)+1, starts[0] + 1 , i+3, ((i+3) - (starts[0]+1))+1)
                    starts = [] #reset list
                    startFound = False 
                    codonFound = True

                if codonFound == False and codon in self.stopsC: #if found stop but not have not found start
                    if ((i+3) - 1)+1 >= self.mGene: #check if codon is appropriate length
                        self.orfSaver((f%3)+1, 1, i+3, ((i+3) - 1)+1) 
                    starts = [] #reset
                    codonFound = True

                #code for reverse sequence begins. Pretty much the same thing as the first except with adjusted start and stop positions to accomodate for the reverse complement.
                if rCodon in self.startsC:
                    rStarts.append(i)
                    rStartFound = True
                if rCodon in self.stopsC and rStartFound == True:
                    #print(head)
                    #print(rStarts)
                    end = len(revSequence) - rStarts[0]
                    length = (len(revSequence) - rStarts[0]) - (len(revSequence) - (i+2)) +1
                    if length >= self.mGene:
                        self.orfSaver(-((f%3)+1), (len(revSequence) - (i+2)), end, length)
                    rStarts = []
                    rStartFound = False
                    rCodonFound = True

                if rCodonFound == False and rCodon in self.stopsC:
                    length = (len(revSequence) - (len(revSequence) - (i - 2))) +1
                    if length >= self.mGene:
                        self.orfSaver(-((f%3)+1), len(revSequence) - (i - 2), len(revSequence), length)
                    rStarts = []
                    rCodonFound = True

        return self.orfDict
        
    def orfSaver(self, frame, start, stop, length):
        '''
        def orfSaver simply saves the orfs and its corresponding properties given to it from orfFinder.
        '''
        self.orfDict.append([frame, start, stop, length])


class NucParams:
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
    dnaNucleotides = { #instantiation of nucleic acid comp
        'A': 0, 'T': 0, 'C':0, 'G':0
    }



    def __init__ (self, fasta):
        '''
        Here is my init method. I personally like to keep all of my methods self contained, so I decided to make this
        accept the fasta sequences here. It's just a way I like doing things I suppose, since this program isn't really 
        designed to take any other input besides genome sequences. 
        '''
        self.aaComp = { #instantiation of aa comp
            value:0 for key, value in NucParams.rnaCodonTable.items()
        }
        self.codonComp = { #instantiation of codon comp
            key:0 for key in NucParams.rnaCodonTable.keys()
        }
        self.nucComp = { #instantiation of nucleic acid comp
            'A': 0, 'T': 0, 'C':0, 'G':0, 'U':0, 'N': 0
        }

        # instantiate fasta reader and iterate through all headers/sequences.
        self.faReader = FastAreader(fasta)
        for head, seq in self.faReader.readFasta():
            self.addSequence(seq)

    def formatSequence(self, sequence):
        '''
        I don't really know if you guys are going to be throwing any oddball sequences so I wrote this to minimize
        the potential of losing points. That way if you guys put in some kind of invalid character it will be removed 
        and the program will still function. If this were an actual program I was creating for research purposes, I would
        definitely delete it due to how much it increases runtime. I'm sure it an be optimized in some way.
        '''
        upperSeq = sequence.upper()
        noWspace = upperSeq.replace(" ", "")
        mySequence = list(noWspace)
        for n, i in enumerate(mySequence):
            if mySequence[n] not in self.nucComp.keys():
                del mySequence[n]
        correctInput = ''.join(mySequence)
        return correctInput
        
    def addSequence (self, seq):
        '''
        Add sequence function. Literally everything is done here. I did it this way to minimize runtime. The way I had it
        set up previously (see code at end) greatly increased runtime since I was basically going to make it run thorough
        the entire library of headers/sequences for each function which is ridiculous. I digress. 
        Limitations of this is that if the sequence is not i%3=0, then it will absolutely not work. But that wasn't
        asked of us in the lab.
        '''

        #this cleans the sequence and returns a nice fresh one with no invalid characters. 
        #greatly increases runtime though.
        nSeq = self.formatSequence(seq)

        #create a new sequence from given sequence that's RNA. Previously I had this set up so it would "check"
        # if the sequence is RNA or DNA, but now everything is RNA. 
        rnaSeq = nSeq.replace("T", "U")
        #create list where each item is a length 3 string (codon)
        rseq = [rnaSeq[i:i+3] for i in range(0, len(rnaSeq), 3)]

        #For each item in the list, check if it matches the key of rnadictionary. If it does, increase value matching its key.
        #  Discard if it doesn't.
        for rCodon in rseq:
            if rCodon in NucParams.rnaCodonTable.keys():
                self.aaComp[NucParams.rnaCodonTable[rCodon]] += 1
                self.codonComp[rCodon] += 1

        #For each item in the sequence, check if its a nucleotide. (there's likely no possible way
        # it won't be due to self.formatSequence())
        for nucleotide in nSeq:
            if nucleotide in self.nucComp.keys():
                self.nucComp[nucleotide] +=1
        
        
    def aaComposition(self):
        return self.aaComp
        
    def nucComposition(self):
        return self.nucComp
        
    def codonComposition(self):
        return self.codonComp

    def nucCount(self):
        return sum(self.nucComp.values())

    #I don't know why I thought I needed this but I kept it in here anyways.
    def codonCount(self):
        return sum(self.codonComp.values())

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class ProteinParam :
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

    def __init__ (self, protein):
        '''
        Initialization of the proteinParams class. Takes any string input from main, formats the input for use by other functions in this program, and determines the amount of each char occurance and stores each char
        as a key, with the amount of occurances as the value.
        '''
        self.protein = self.formatInput(protein)
        self.aaComp = {}
        self.aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        #Counts all occurances of each char item in self.aa in the formatted input string
        for self.aa in self.aa:
            if self.aa in self.protein: #if char in inputstring matches char in list
                num = self.protein.count(self.aa) #count number of instances
                self.aaComp.update({self.aa:num}) #update dictionary with new value
            else: # else if char not in input string
                self.aaComp.update({self.aa:0}) #update dictionary with zero
        self.Cystine = True

    def formatInput (self, sequence):
        '''
        Formats the input string for use by this program. Makes all letters uppercase, replaces all whitespace, and deletes all chars that are not keys in aa2mw.
        '''
        upperSeq = sequence.upper()
        noWspace = upperSeq.replace(" ","")
        myProtein = list(noWspace)
        for n, i in enumerate(myProtein):
            if myProtein[n] not in ProteinParam.aa2mw.keys():
                del myProtein[n]
        correctInput = ''.join(myProtein)
        return correctInput

    
    def aaCount (self):
        '''
        Returns the length of the formatted string.
        '''
        return len(self.protein)



    def pI (self):  
        '''
        Binary search implementation of pI. Since midCharge will never actually be zero and only approach zero, the number midCharge must be rounded. Increasing the y in the round(x, y)
        function will give a more accurate result, but the runtime of this algorithm is mainly dependent on y and will increase the higher that number is (and I mean that by effectively increasing the "length" 
        of the "array" which midCharge is disguised as- the true runtime of this is still log(n).)
        ''' 
        L = 0.0
        R = 14.01
        pI = 0
        midpH = 0

        while L < R:
            midpH = (L+R)/2 #determine arbitrary midpoint
            midCharge = self._charge_(midpH) # find charge at that midpoint
            if round(midCharge, 2) == 0.00: # if midCharge = target return midpH
                pI = midpH
                return pI
            elif round(midCharge, 2) > 0.00: # else throw out left half
                L = midpH
            elif round(midCharge, 2) < 0.00: # else throw out right half
                R = midpH

    def aaComposition (self):
        '''
        See __init__ for composition of the dictionary. 
        '''
        return self.aaComp

    def _charge_ (self, pH):
        '''
        Function for determining charge of a protein using the values of aa2chargePos and aa2chargeNeg
         and matching keys between aa2chargePos/aa2chargeNeg and aaComp, and the values of aaComp. Also takes into account N terminus and C terminus. 
        '''
        pnetCharge = (10**ProteinParam.aaNterm) / (10**ProteinParam.aaNterm + 10**pH) #instantiate pnet with value for N term
        nnetCharge = (10**pH) / (10**ProteinParam.aaCterm + 10**pH) #instantiate nnet with value for C term
        
        for self.aa in ProteinParam.aa2chargePos: #iterate through every aa
            n = self.aaComp[self.aa] # get value for key (number of occurance)
            pnetCharge += n * ((10**ProteinParam.aa2chargePos[self.aa]) / ( (10**ProteinParam.aa2chargePos[self.aa]) + 10**pH ))
        for self.aa in ProteinParam.aa2chargeNeg:
            n = self.aaComp[self.aa] #get value for key (number of occurance)
            nnetCharge += n * ((10**pH) / ( (10**ProteinParam.aa2chargeNeg[self.aa]) + 10**pH ))

        netCharge = pnetCharge - nnetCharge
        
        return netCharge

    def molarExtinction (self):
        '''
        Determines amount of light a protein absorbs at 280nm by multiplying the values of aaComp and aa2abs280
        so long as their keys match. Provided additional functionality for "reducing conditions" which can be adjusted in __init__ 
        under self.Cystine.
        '''
        E = 0
        if self.Cystine == True: 
            for self.aa in ProteinParam.aa2abs280:
                E += self.aaComp[self.aa] * ProteinParam.aa2abs280[self.aa]
        else: # if self.Cystine = false
            for self.aa in ProteinParam.aa2abs280:
                if self.aa == "C": # if self.aa = Cystine
                    E += 0 # skip
                else: 
                    E += self.aaComp[self.aa] * ProteinParam.aa2abs280[self.aa] #else continue normally
        return E
        

    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        Method for determining molecular weight of a protein. Function simply iterates through each char in the formatted input string,
        and if the char matches a key in aa2mw, will add the value in aa2mw subtracted by the value of mwH2O.
        '''
        tMW = ProteinParam.mwH2O
        for n, i in enumerate(self.protein): #iterate through each char within input string. There are likely more pythonic ways to do this, but I decided to go with this.
            if self.protein[n] in ProteinParam.aa2mw.keys():    
                tMW = tMW + (ProteinParam.aa2mw[self.protein[n]] - ProteinParam.mwH2O)
        return tMW
