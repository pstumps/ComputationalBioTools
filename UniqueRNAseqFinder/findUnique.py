#!/usr/bin/env python3
# Name: Patrick Stumps (pstumps)
# Group Members: None

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sequenceAnalysis import FastAreader

'''
findUnique.py
Here is my implementation of findUnique.py. In this implementation, I create a tRNA object that is composed of a header, a sequence, a powerset 
for its own respective sequence, a unique set for its own respective sequence composed of the difference between the union of all powersets minus 
its own powerset, an essential set composed of all of the unique and essential elements for its respective set, and lastly a dictionary of all
the unique and essential elements where the key is the elements index in its respective sequence and the value is the element. 

Usage:
input: python3 findUnique.py [inputfasta]
output: unique and essential elements
'''

class findUnique:
    '''
    class findUnique
    Here is my findUnique class. It is here that the program instantiates a list where it will save all of the rna objects (header, sequence, powerset, etc.)
    Before this program can do anything, it must create a powerset of all of the sequences in the fastAfile. It must do so in order to determine the unique
    set for the rnas sequence.
    '''

    rnaSetlist = list()

    def __init__(self, header, sequence):
        '''
        __init__
        The init method is where an rna object instantiates its parameters. Once these parameters have been established, 
        it will append itself to the list of rnas. The dictionary of unique and essential elements is probably not necessary, but I believe it would be useful
        if someone wanted to "reach" into a sequence and find a specific element at a specific index instead of popping them out of the essential set.
        This could probably be done using pythons built in find() function but the dictionary saves them the trouble.
        '''
        self.header = header.replace(" ", "")
        s1 = sequence.replace("_", "")
        self.sequence = s1.replace("-", "") #remove invalid characters from alignment
        self.powerSet = set()
        self.unique = set()
        self.eSet = set()
        self.eDictionary = {}
        for i in range(len(self.sequence)): 
            for j in range(i+1, len(self.sequence)+1):
                self.powerSet.add(self.sequence[i:j]) 
        self.rnaSetlist.append(self) # append self object consisting of all of the above parameters to the list of objects
    
    def uniqueFinder(self):
        '''
        uniqueFinder
        Here's where the magic happens. First the program creates a set of nonUnique elements by creating the union of all powersets except the objects own
        powerset. It then determines the elements unique to the objects own sequence by computing the difference of its own powerset minus the union of all
        powersets. The program then determines if each unique element is essential by creating a slice of the unique element from the beginning to the length
        of the unique element and iterating through, increasing the start of the slice which therefore decreases the length of the slice, then adding the slice
        to a set (I believe this is essentially the extension method we talked about in class). Finally, the program creates the dictionary of essential elements
        (discussed earlier).
        '''
        nonUnique = set()

        #get unique elements
        for rna in self.rnaSetlist:# for objects in object list
            if rna.powerSet is not self.powerSet: #if object powerset is not same as current object powerset
                nonUnique = nonUnique.union(rna.powerSet) #create union of powerset

        self.unique = self.powerSet.difference(nonUnique) #create unique set from difference between self object power set and nonunique set
        #print(self.unique)
        #get essential elements

        for uniqueElement in sorted(self.unique, key=len): #for element(string) in unique set
            uLen = len(uniqueElement) # instantiate length
            if not any(uniqueElement[i:i+j+1] in self.eSet for i in range(uLen) for j in range(uLen - i)): #if slice from i to element length - i is not in set
                self.eSet.add(uniqueElement) # add element
        
        #make essential element dictionary
        for k in self.eSet:
            index = self.sequence.find(k)
            self.eDictionary[index] = k

        return self.eDictionary
    '''
    def getPowerSets(self, sequence):
    '''
        #getPowerSets
        #I originally intended for this to be a generator in my init method but python didn't like that as it limited the amount of functions I could perform
        #(len(), string[i:j] - not subscriptable) that are necessary for my program to work. So I went with this.
    '''
        for i in range(len(self.sequence)): 
            for j in range(i+1, len(self.sequence)+1):
                self.powerSet.add(self.sequence[i:j])
    '''

def main():
    '''
    main
    Execute all functions in order for proper output to stdout.
    '''
    
    fReader = FastAreader(sys.argv[1])
    for head, seq in fReader.readFasta():
        print("Generating powerset for: {}".format(head))
        findUnique(head, seq)
        #fa.getPowerSets(seq)
    
    for rna in sorted(findUnique.rnaSetlist, key = lambda rna:rna.header):
        uniqueAndEssential = rna.uniqueFinder()
        print(rna.header)
        print(rna.sequence)
        for index, seq in sorted(uniqueAndEssential.items(), key = lambda x:x[0]):
            print('.'*index, seq, sep ='')

main()
