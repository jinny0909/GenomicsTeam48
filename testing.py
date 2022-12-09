# Computational Genomics Final Project
# Team 48
"""This file is for testing purposes.  It can be run to replicate the results in the writeup, and it uses all
 the methods created for this project.  The three tests are for exact protein alignment, inexact protein alignment
 and checking that terminators do not match."""

from ctypes import alignment
from bwa import *
from proteinToDNA import *
from proteinIndex import *
from openRF import *
from blosum import *
import random 

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

def test(p_name, protein, mismatch):
    print("Protein Name: ")
    print(p_name)
    DNA = proteinToDNA(protein)
    tests = getRF(DNA)
    passed_test, passed_protein = [], []

    for test1 in tests:
        output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs, mismatch)
    for o in output:
        passed_test.append(test1)
        passed_protein.append(proteins[protein_id[o]])
        print("Protein match: %s" % protein_id[o])

    if len(passed_protein) > 1 and len(passed_test) > 1: 
        print("total passed")
        print(len(passed_protein))
        protein_index = FMIndex(passed_protein[0] + '$', 10, 5)
        output = protein_index.query(passed_test[0])
        print(output)
        assert(len(output) == 1)
        assert(output != (-1,-1))

        perfect_score, alignment_score = localAlignment(passed_protein[0], passed_test[0], -1, -1,  True)
        print(perfect_score)
        print(alignment_score)
        if mismatch == 0: 
            assert(perfect_score == alignment_score)

# Exact Match Tests 
name_list = [' 3-isopropylmalate dehydrogenase ', ' 2-iminoacetate synthase ', ' 1,4-dihydroxy-2-naphthoyl-CoA synthase ', 
' 2,5-diketo-D-gluconic acid reductase B ', ' 2-hydroxy-3-oxopropionate reductase ' ]
for i in range(5):
    p_name = name_list[i]
    protein_seq = proteins[p_name]
    protein_seq = protein_seq[5:40]
    print("Test %d" %(i+1))
    test(p_name, protein_seq, 0)
    




