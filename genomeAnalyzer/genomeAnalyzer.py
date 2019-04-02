#!/usr/bin/env python3
#Name: Patrick Stumps(pstumps)


import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sequenceAnalysis import NucParams

def main ():
    '''
    main function. Constructs libraries for aa composition, nucleotide composition and codon composition.
    '''

    myNuc = NucParams(sys.argv[1])

    aacDict = myNuc.aaComposition() #returns dictionary of form aa: #ofaaoccurance
    naDict = myNuc.nucComposition() #returns dictionary of form base: #ofbase
    ccDict = myNuc.codonComposition() #returns dictionary of form codon: #ofcodonoccourance

    #Determine sequence length
    print("sequence length = {:.2f} Mb".format(myNuc.nucCount()*0.000001))
    print()

    #sequence GC content
    print("GC content = {:.1f}%".format((naDict['G'] + naDict['C'])*100/myNuc.nucCount()))
    print()
    
    #for all keys/values in rna dictionary, iterate through sorted keys then values
    for nuc, aa in sorted(NucParams.rnaCodonTable.items(), key=lambda x: x[1]+x[0]):
        try:
            val = ccDict[nuc] / aacDict[NucParams.rnaCodonTable[nuc]] #divide the amount of nucleotide by amount of amino acid
        except ZeroDivisionError:
            val = 0 # if aadict[key] is zero return zero, negates divbyzero error
        except KeyError:
            val = 0 # if aadict[key] is zero return zero
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, ccDict[nuc]))
            #print function in order : nucleotide, amino acid (rnacodontable key), codon frequency, codon count.
if __name__ == "__main__":
    main()