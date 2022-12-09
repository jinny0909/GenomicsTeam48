# Computational Genomics Final Project
# Team 48
"""This file is for testing purposes.  It can be run to replicate the results in the writeup, and it uses all
 the methods created for this project.  The three tests are for exact protein alignment, inexact protein alignment
 and checking that terminators do not match."""

from bwa import *
from proteinToDNA import *
from proteinIndex import *
from openRF import *
from blosum import *

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

#NA = DNA[0:2] + 'A' + DNA[3:] 
print(DNA)

# open reading frame
passed_test = []
passed_protein = []
tests = getRF(DNA)
for test1 in tests:
    output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs, 1)
    for o in output:
        passed_test.append(test1)
        passed_protein.append(proteins[protein_id[o]])
        print("Protein match: %s" % protein_id[o])

print(passed_test)
print(passed_protein)
# verify using FM
protein_index = FMIndex(passed_protein[0] + '$', 10, 5)
output = protein_index.query(passed_test[0])
print(output)

print(len(passed_test[0]))
localAlignment(passed_protein[0], passed_test[0], -1, -1,  True)
