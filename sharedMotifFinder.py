#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given a collection of DNA string in fasta format program finds the largest
common motif among all strings. Returns this common motif 

Created on Wed Jul 15 10:41:20 2020
@author: bram
"""

# Import modules needed
import sys

# Function comapres two sequences and finds common motifs between them
def compare2(string1: str, string2: str):
    """
    Given two DNA sequences returns a list containing all the regions
    longer than two bases that are the same between the two sequences
    """
    same = set() # List to hold motifs
    
    # Variables to hold longer and shorter sequences and their lengths
    longSeq = ''
    shorter = 0 
    shortSeq = ''
    
    # Determine lengths of strings relative to each other. 
    #   Doesn't matter if they are the same length or not
    if len(string1) < len(string2):
        shorter = len(string1)
        shortSeq = string1
        longSeq = string2
    else:
        longSeq = string1
        shortSeq = string2
        shorter = len(string2)
    
    for i in range(2, shorter + 1):
        for j in range (0, shorter):
            if shortSeq[j:j+i] in longSeq:
                same.add(shortSeq[j:j+i])
    return [motif for motif in same if len(motif) > 1]

def commonMotif(motifList, string1: str):
    """
    Given a list of motifs and a DNA sequence returns a motif list 
    containing those motifs, out of the motifs in the list, that were
    found in the sequence
    """
    common = [] # List to hold motifs
    for motif in motifList:
        if motif in string1:
            common.append(motif) # If motif is present in sequence add to list
    return common

# List to hold IDs and sequences from raw fasta file 
DNAStringArray = []

# Open and parse input file
with open(sys.argv[1]) as f:
    lines = [item.strip() for item in f.readlines() if not item == '']
    tempList = [] # list for sequences associated with each ID. Solves sequences on multiple lines issue.
    for line in lines:
        # All Lines with sequence IDs start with '>'
        if line.startswith('>'):
            # If tempList contains sequences when new ID reached add 
            #   ID and sequences in tempList to dictionary. Empty tempList
            if not len(tempList) == 0:
                DNAStringArray.append(tempList[1:]) # Add sequences to list
                tempList = [line,] # Empty List, Add new ID to it 
            else: # Add seq ID to list if it's the first sequence ID
                tempList = [line,]
        # Add sequences to the tempList if no ID is encountered
        else:
            tempList.append(line)
    # Handles sequence for last ID
    else:
        DNAStringArray.append(tempList[1:]) # Add sequences to list

# Save first two sequences from 
seq1, seq2, *rest = [item for item in DNAStringArray]
seq1 = ''.join(seq1) 
seq2 = ''.join(seq2)

# Find common motifs among the first two sequences to start the comparison
motifs = compare2(seq1, seq2)

# Filter based on similarities in the other sequences 
for i in rest:
    seq = ''.join(i)
    motifs = commonMotif(motifs, seq)

print(motifs) # Motifs common among all sequences 
print(max(motifs, key=len)) # Return the longest item