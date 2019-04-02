#!/usr/bin/env python3
# Name: Patrick Stumps (pstumps)
# Group Members: None
'''
findORFs.py
Here is my implementation of findORFs.py. To be honest the main algorithm for the program isn't really here but this is where the magic starts. 
It should really be called something like "main.py" but I digress. This program takes all the commandline input and passes it to the
the findORF class for instantiation. It will also create the correct output file and print all of the information given to it by findorf in the appropriate
order given.

Input: fa file, commandline parameters
Output: fa outfile with all orfs
'''
######################################################################
##
# CommandLine
######################################################################
##
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
from sequenceAnalysis import orfFinder, FastAreader

class CommandLine(): 
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various
    argument options, a standard usage and help.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    def __init__(self, inOpts = None) : 
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        self.parser = argparse.ArgumentParser(description = 'Program prologe - a brief description of what this thing does', 
                                            epilog = 'Program epilog- some other stuff you feel compelled to say', #default is True
                                            add_help = True, #default prefix_chars = '-',
                                            usage = '%(prog)s[options] -option1[default] <input>output')
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name')
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (0, 5, 10, 100,200,300,500,1000), default=100, action = 'store', help='minimumGene length- choices are 100, 200, 300, 500, 1000')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'], nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows mult iple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)
        
def main(inCL=None): 

    #Find some genes.
    '''
    Main function. Initializes all items from args and passes them to findORFs in sequenceAnalysis.py. This will also type out all the sorted correct output
    into the outfile given in the namespace.
    '''
    if inCL is None:
        myCommandLine = CommandLine() 
    else :
        myCommandLine = CommandLine(inCL)
    print(myCommandLine.args) 
    inFile = myCommandLine.args.inFile
    outFile = myCommandLine.args.outFile
    longestGene = myCommandLine.args.longestGene
    minGene = myCommandLine.args.minGene
    startCodons = myCommandLine.args.start
    stopCodons = myCommandLine.args.stop

    fReader = FastAreader(inFile)
    orfOut = open(outFile, 'w')
    for head, seq in fReader.readFasta():
        Orfer = orfFinder(longestGene, minGene, startCodons, stopCodons)
        orfOut.write(head+ "\n")
        Orfer.findOrf(head, seq)
        for frame, start, stop, length in sorted(Orfer.orfDict, key =lambda kv: (-kv[3], -kv[1])): # ensures everything will be printed out in descending order from gLength values then keys.
            orfOut.write("{:+d} {:>5d}..{:>5d} {:>5d} \n".format(frame, start, stop, length))
    #######
#if __name__ == "__main__": main(['tass2.fa', 'tass2ORFdata-ATG-100.txt', '--longestGene']) # delete the list when you want to run STDIN
if __name__ == "__main__": main()
'''
testing purposes
def main():
    myReader = FastAreader('tass2.fa')
    for head, seq in myReader.readFasta():
        orf = orfFinder()
        orf.findOrf(seq)
'''