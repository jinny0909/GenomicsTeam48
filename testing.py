# Computational Genomics Final Project
# Team 48
"""This file is for testing purposes.  It can be run to replicate the results in the writeup, and it uses all
 the methods created for this project.  The three tests are for exact protein alignment, inexact protein alignment
 and checking that terminators do not match."""

from bwa import *
from proteinToDNA import *
from proteinIndex import *
from openRF import *

# dataset - 100 protein sequences for E. Coli chosen by name (alphabetical order)
dataset = 'uniprot-compressed_true_download_true_format_fasta_query_accession_3-2022.12.08-22.52.25.99.fasta'
protein_id = dict()
proteins = dict()
count = 1
reference = ""
with open(dataset, 'r') as file:
    temp = file.readline()
    while temp != '':
        if temp[0:2] == '>s':
            name = temp.split("ECOLI", 2)
            name = name[1].split("OS", 2)
            name = name[0]
            temp = file.readline()
            protein_id[count] = name
            count += 1
        else:
            sequence = ''
            while temp[0:2] != '>s' and temp != '':
                sequence = sequence + temp.rstrip()
                temp = file.readline()
            proteins[name] = sequence
            reference = reference + sequence + '&'

SA, C, O, O_rev, n, alphabet, proteinIDs = precalculation(reference + '$')

# Part 1 - Exact matching
"""Test 1"""
# protein name: 23S rRNA 5-hydroxycytidine C2501 synthase
protein_sequence = "QIVLARELNLDQIRAIHQATDATIEF"
# change to DNA
DNA = proteinToDNA(protein_sequence)

# open reading frame
tests = getRF(DNA)
for test1 in tests:
    output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs)
    for o in output:
        print("Protein match: %s" % protein_id[o])

# verify using FM
protein_index = FMIndex(reference + '$', 20, 5)
output = protein_index.query(test1)
print(output)


