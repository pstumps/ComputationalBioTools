#!/usr/bin/env python3
#Name: Patrick Stumps(pstumps)
#Group Members: none
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fileinput
import argparse
import numpy as np

from sequenceAnalysis import NucParams

def main ():
    '''
    Same thing as my main function in GenomeAnalyzer.py except it takes two arguements from stdin, determines sequence length for both, GC content for both,
    then takes the log2 ratio of both genomes 1 and two. The output of the last lines is also explained in the print function.
    '''

    genome1 = NucParams(sys.argv[1])
    genome2 = NucParams(sys.argv[2])

    aacDict1 = genome1.aaComposition() #returns dictionary of form aa: #ofaaoccurance
    naDict1 = genome1.nucComposition() #returns dictionary of form base: #ofbase
    ccDict1 = genome1.codonComposition() #returns dictionary of form codon: #ofcodonoccourance

    aacDict2 = genome2.aaComposition()
    naDict2 = genome2.nucComposition()
    ccDict2 = genome2.codonComposition()

    print("Genome 1 sequence length = {:.2f} Mb".format(genome1.nucCount()*0.000001))
    print()
    print("Genome 2 sequence length = {:.2f} Mb".format(genome2.nucCount()*0.000001))
    print()

    print("Genome 1 GC content = {:.1f}%".format((naDict1['G'] + naDict1['C'])*100/genome1.nucCount()))
    print()
    print("Genome 2 GC content = {:.1f}%".format((naDict2['G'] + naDict2['C'])*100/genome2.nucCount()))
    print()

    g1GC = (naDict1['G'] + naDict1['C'])*100/genome1.nucCount()
    g2GC = (naDict2['G'] + naDict2['C'])*100/genome2.nucCount()

    print("GC content log2 ratio of Gnome 1 and Genome 2 = {:.1f}%".format(np.log2(g1GC/g2GC)))
    print()
    
    for nuc, aa in sorted(NucParams.rnaCodonTable.items(), key=lambda x: x[1]+x[0]):
        try:
            val = ccDict1[nuc] / aacDict1[NucParams.rnaCodonTable[nuc]]
            val2 = ccDict2[nuc] / aacDict2[NucParams.rnaCodonTable[nuc]]
            
        except ZeroDivisionError:
            val = 0
            val2 = 0
        except KeyError:
            val = 0
            val2 = 0

        try:
            log2ratio = np.log2(val/val2)
        except ZeroDivisionError:
            log2ratio = 0
        except KeyError:
            log2ratio = 0
        
        print ('{:s} : {:s} | {:5.1f} {:5.1f} | ({:6d} {:6d}) | {:6f}'.format(nuc, aa, val*100, val2*100, ccDict1[nuc], ccDict2[nuc], log2ratio))
    print('AA : Codon | Genome 1 codon frequency, Genome 2 codon frequency | (Genome 1 codon count, Genome 2 codon count) | log2 ratio of codon frequencies')
if __name__ == "__main__":
    main()