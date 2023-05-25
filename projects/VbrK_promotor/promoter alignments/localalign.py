"""Finds clusters of best alignments of 2 given sequences"""

import re
import os
import math
import sys

class MinHeap:
  
    def __init__(self, maxsize):
        self.maxsize = maxsize
        self.size = 0
        self.Heap = [0]*(self.maxsize + 2)
        self.Heap[0] = (-100, -100)
        self.FRONT = 1

    def parent(self, pos):
        return pos//2

    def leftChild(self, pos):
        return 2 * pos

    def rightChild(self, pos):
        return (2 * pos) + 1

    def isLeaf(self, pos):
        return pos*2 > self.size

    def swap(self, fpos, spos):
        self.Heap[fpos], self.Heap[spos] = self.Heap[spos], self.Heap[fpos]
  
    # Function to heapify the node at pos
    def minHeapify(self, pos):
  
        # If the node is a non-leaf node and greater
        # than any of its child
        if not self.isLeaf(pos):
            if (self.Heap[pos][0] > self.Heap[self.leftChild(pos)][0] or 
               self.Heap[pos][0] > self.Heap[self.rightChild(pos)][0]):
  
                # Swap with the left child and heapify
                # the left child
                if self.Heap[self.leftChild(pos)][0] < self.Heap[self.rightChild(pos)][0]:
                    self.swap(pos, self.leftChild(pos))
                    self.minHeapify(self.leftChild(pos))
  
                # Swap with the right child and heapify
                # the right child
                else:
                    self.swap(pos, self.rightChild(pos))
                    self.minHeapify(self.rightChild(pos))
  
    # Function to insert a node into the heap
    def insert(self, element):
        self.size+= 1
        self.Heap[self.size] = element
  
        current = self.size
  
        while self.Heap[current][0] < self.Heap[self.parent(current)][0]:
            self.swap(current, self.parent(current))
            current = self.parent(current)
        
        if self.size > self.maxsize:
            self.remove()

    # Function to print the contents of the heap
    def Print(self):
        for i in range(1, (self.size//2)+1):
            print(" PARENT : "+ str(self.Heap[i])+" LEFT CHILD : "+ 
                                str(self.Heap[2 * i])+" RIGHT CHILD : "+
                                str(self.Heap[2 * i + 1]))
  
    # Function to build the min heap using
    # the minHeapify function
    def minHeap(self):
  
        for pos in range(self.size//2, 0, -1):
            self.minHeapify(pos)
  
    # Function to remove and return the minimum
    # element from the heap
    def remove(self):
  
        popped = self.Heap[self.FRONT]
        self.Heap[self.FRONT] = self.Heap[self.size]
        self.size-= 1
        self.minHeapify(self.FRONT)
        return popped

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

class Matrix:
    def __init__(self, row, col):
        self.row = row
        self.col = col
        self.row_head = list(range(1, row+1))
        self.col_head = list(range(1, col+1))
        self.matrix = []
        for i in range(row):
            self.matrix.append([])
        for i in range(row):
            for j in range(col):
                self.matrix[i].append(0)
    
    def max(self):
        self.max_entry = 0
        for i in range(self.row):
            for j in range(self.col):
                if i == j:
                    continue
                if self.matrix[i][j] > self.max_entry:
                    self.max_entry = self.matrix[i][j]
        return self.max_entry
        
    def min(self):
        self.min_entry = 100
        for i in range(self.row):
            for j in range(self.col):
                if i == j:
                    continue
                if self.matrix[i][j] < self.min_entry:
                    self.min_entry = self.matrix[i][j]
        return self.min_entry
    
    def display(self):
        cell_width = (os.get_terminal_size().columns // (self.col+1)) - 2
        print("|", "".center(cell_width)[:cell_width], end="|", sep="")
        for i in self.col_head:
            print(str(i).center(cell_width)[:cell_width], end="|")
        else:
            print("")
        for i in range(self.row):
            print("|", str(self.row_head[i]).center(cell_width)[:cell_width], end="|", sep="")
            for j in self.matrix[i]:
                print(str(j).center(cell_width)[:cell_width], end="|")
            else:
                print("")

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

def align_local(seq1, seq2):
    output = 1000
    match = 5
    mismatch = -3
    gap = -3
    
    if seq1.alphabet != seq2.alphabet:
        raise Exception("Alphabet of sequences does not match")
    
    #create matrix
    align_mat = Matrix(len(seq1.seq)+1, len(seq2.seq)+1)
    align_mat.row_head = list(" " + seq1.seq)
    align_mat.col_head = list(" " + seq2.seq)
    
    #initialize 1st row/collumn
    for i in range(align_mat.row):
        align_mat.matrix[i][0] = (0, "e")
    for i in range(align_mat.col):
        align_mat.matrix[0][i] = (0, "e")
    
    #calculate values
    for i in range(1, align_mat.row):
        for j in range(1, align_mat.col):
            if align_mat.row_head[i] == align_mat.col_head[j]:
                align_mat.matrix[i][j] = (align_mat.matrix[i-1][j-1][0] + match, "d")
            else:
                score_dict = {}
                score_dict[align_mat.matrix[i][j-1][0] + gap] = "l"
                score_dict[align_mat.matrix[i-1][j][0] + gap] = "u"
                score_dict[align_mat.matrix[i-1][j-1][0] + mismatch] = "d"
                score_dict[0] = "e"
                
                score = max(score_dict.keys())
                align_mat.matrix[i][j] = (score, score_dict[score])
    
    #find highest values
    min_heap = MinHeap(output)
    for i in range(align_mat.row):
        for j in range(align_mat.col):
            min_heap.insert((align_mat.matrix[i][j][0], i, j))
    
    max_val = []
    for i in range(min_heap.maxsize):
        max_val.append(min_heap.remove())

    #backtrace
    align_sequences = []
    for value in max_val:
        end = value[1]
        i = value[1]
        j = value[2]
        align_seq1 = ""
        align_seq2 = ""
        align_score = value[0]
        
        while align_mat.matrix[i][j][1] != "e":
            if align_mat.matrix[i][j][1] == "d":
                align_seq1 = align_mat.row_head[i] + align_seq1
                align_seq2 = align_mat.col_head[j] + align_seq2
                i -= 1
                j -= 1
            elif align_mat.matrix[i][j][1] == "l":
                align_seq1 = "-" + align_seq1
                align_seq2 = align_mat.col_head[j] + align_seq2
                j -= 1
            elif align_mat.matrix[i][j][1] == "u":
                align_seq1 = align_mat.row_head[i] + align_seq1
                align_seq2 = "-" + align_seq2
                i -= 1
        
        align_sequences.append((i, end, align_score, align_seq1, align_seq2))
    
    return align_sequences

def cluster_align(sequences):
    clusters = {}
    tolerance = 20
    
    for seq in reversed(sequences):
        for best in clusters.keys():
            if seq[0] > (best[0] - tolerance) and seq[1] < (best[1] + tolerance):
                clusters[best][0] += 1
                clusters[best][1] += seq[2]
                break
        else:
            clusters[seq] = [1, seq[2]]
    
    return clusters

def main():
    if len(sys.argv) != 3:
        print("Incorrect arguments given (sequence1, sequence2)")
        quit()
    _, fasta1, fasta2 = sys.argv
    
    seq1 = read_fasta(fasta1)[0]
    seq2 = read_fasta(fasta2)[0]
    
    align = align_local(seq1, seq2)
    cluster = cluster_align(align)
    
    print("Optimal local alignments")
    for i in cluster.keys():
        print("\nTop Score: " + str(i[2]) + " Start: " + str(i[0]) + " End: " + str(i[1]))
        print(str(i[3]) + "\n" + str(i[4]))
        print("with " + str(cluster[i][0]) + " similar alignments and avg score of " + str(cluster[i][1]/cluster[i][0]))

if __name__ == "__main__":
    main()