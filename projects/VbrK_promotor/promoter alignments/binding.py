"""Finds clusters of best alignments of 2 given sequences"""

import re
import os
import math
import sys

class Sequence:
    def __init__(self, seq, header):
        self.id = 0
        self.seq = seq
        self.header = header
        #isolates the genename in the header through regex search
        self.genename = re.findall("\[gene=[^\s]*", header)
        self.alphabet = self.identifyAlphabet()
        self.frequencies = self.identifyFrequencies()

    def identifyAlphabet(self):
        """returns if the fasta file describes a gene or protein"""
        for i in self.seq:
            if i not in ("A", "T", "G", "C", "-"):
                return "protein"
        else:
            return "gene"

    def identifyFrequencies(self):
        """returns dict with absolute and relative frequencies of the alphabet in seq"""
        freq = {}
        length = len(self.seq)
        
        #count absolute frequency of letters
        for i in self.seq:
            if i in freq:
                freq[i] += 1
            else:
                freq[i] = 1
        
        #calculate relative frequency
        for i in freq.keys():
            freq[i] = (freq[i], freq[i]/length)
        
        return freq

def read_fasta(path):
    """reads the fasta file at path and returns a fasta object based of said file"""
    new_sequences = []
    
    with open(path, "r") as file:
        raw_data = file.read()
        file.close()
    
    #splits the data in case of multifasta
    raw_data = raw_data.split("\n")
    header = ""
    sequence = ""
    for i in raw_data:
        if i and i[0] == ">":
            if header:
                new_sequences.append(Sequence(sequence, header))
            header = i
            sequence = ""
        else:
            sequence += i
    else:
        new_sequences.append(Sequence(sequence, header))
    
    return new_sequences

def mismatches(seq1, seq2):
    """counts number of mismatches in two equal length sequences (gaps dont count)"""
    
    #check if sequences are equal length
    if len(seq1) != len(seq2):
        return -1
    
    mismatches = 0
    
    #check each position in string for a mismatch
    for i in range(len(seq1)):
        #ignore gaps
        if seq1[i] == "_" or seq2[i] == "_":
            continue
        
        if seq1[i] != seq2[i]:
            mismatches += 1
    
    return mismatches

def seq_search(sequence, gap_min, gap_max, s1, s2, allowed_mismatch):
    """searches sequence for all ocurrences of s1 and s2"""
    
    #list of matching substrings and the matching query
    return_seqs = []
    
    for gap_num in range(gap_min, gap_max+1):
        #construct complete query with desired gapcount
        query = s1 + "_"*gap_num + s2
        
        #check all possible substrings of sequence for query
        for i in range(len(sequence) - len(query) + 1):
            subject = sequence[i:i+len(query)]
            mismatch = mismatches(subject, query)
            if mismatch <= allowed_mismatch:
                return_seqs.append((subject, query, i, mismatch))
    
    return return_seqs

def main():
    if len(sys.argv) != 2:
        print("Incorrect arguments given (sequence)")
        quit()
    _, fasta = sys.argv
    
    #define parameters
    seq = read_fasta(fasta)[0]
    s1 = "TTCTAAT"
    s2 = "TTCATCG"
    #amount of gaps between s1 and s2
    min_gap = 1
    max_gap = 6

    #allowed mismatches of s1+s2 from sequence excluding gaps
    allowed_mismatch = 5
    
    #search matches
    seqs = seq_search(seq.seq, min_gap, max_gap, s1, s2, allowed_mismatch)
    
    #output
    print("Found matches " + str(len(seqs)) + " of " + s1 + " and " + s2 + ":")
    for i in seqs:
        print("Index:" + str(i[2]) + "," + str(i[2] + len(i[1])) + " Mismatch: " + str(i[3]))
        print(i[0])
        print(i[1])
        print()

if __name__ == "__main__":
    main()