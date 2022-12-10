# Computational Genomics Final Project
# Team 48
"""This file is for testing purposes.  It can be run to replicate the results in the writeup, and it uses all
 the methods created for this project.  The three tests are for exact protein alignment, inexact protein alignment
 and checking that terminators do not match.  It is basically formatted in the way of our full tool"""
import time
from ctypes import alignment
from bwa import *
from proteinToDNA import *
from proteinIndex import *
from openRF import *
from blosum import *
import random
from parse_data import *

# dataset - 100 protein sequences for E. Coli chosen by name (alphabetical order)
dataset = 'uniprot-compressed_true_download_true_format_fasta_query_accession_3-2022.12.08-22.52.25.99.fasta'

# parse data
reference, protein_id, proteins = parse_data(dataset)
SA, C, O, O_rev, n, alphabet, proteinIDs = precalculation(reference + '$')
protein_index = FMIndex(reference + '$', 10, 5)


def test(DNA, mismatch, check):
    tests = getRF(DNA)
    passed_test, passed_protein = [], []

    for test1 in tests:
        start = time.time()
        output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs, mismatch)
        end = time.time()
        # print('BWA elapsed time (millisec):')
        # print((end-start)*1000)
        for o in output:
            passed_test.append(test1)
            passed_protein.append(proteins[protein_id[o]])
            print("Protein match: %s" % protein_id[o])

            dp_start = time.time()
            perfect_score, alignment_score = localAlignment(passed_protein[0], passed_test[0], -1, -1,  True)
            dp_end = time.time()
            print("Perfect Alignment Score:")
            print(perfect_score)

            print('Alignment Score')
            print(alignment_score)
            print('DP elapsed time (millisec):')
            print((dp_end-dp_start)*1000)
            print("\n")


    if check == 0:
        if len(passed_protein) > 1 and len(passed_test) > 1:
            print("total passed")
            print(len(passed_protein))
            output = protein_index.query(passed_test[0])
            print(output)
            assert(len(output) == 1)
            assert(output != (-1, -1))

      
            if mismatch == 0:
                assert(perfect_score == alignment_score)

name_list = [' 3-isopropylmalate dehydrogenase ', ' 2-iminoacetate synthase ', ' 1,4-dihydroxy-2-naphthoyl-CoA synthase ', 
' 2,5-diketo-D-gluconic acid reductase B ', ' 2-hydroxy-3-oxopropionate reductase ' ]

dnas = []
# Tests for Exact Matching
print("Test Exact")
print("1 MATCH SHOULD OCCUR FOR EACH TEST")
print("-----------------------------------------")
for i in range(5):
    p_name = name_list[i]
    protein_seq = proteins[p_name]
    protein_seq = protein_seq[5:40]
    print("Test %d" % (i+1))
    print("Protein Name: ")
    print(p_name)
    DNA = proteinToDNA(protein_seq)
    DNA = transcribe(DNA)
    dnas.append(DNA)
    test(DNA, 0, 0)

# Check that mismatch gives no Exact Matches
print("\nNO MATCHES SHOULD OCCUR FOR FOLLOWING MISMATCHES")
print("-----------------------------------------")
print("Test 1")
test(dnas[0][0:10] + 'C' + dnas[0][10:], 0, 1) # insertion
print("Test 2")
test(dnas[1][0:9] + dnas[1][10:], 0, 1) # deletions
print("Test 3")
test(dnas[2][0:9] + 'A' + dnas[2][10:], 0, 1) # substitution
print("Test 4")
test(dnas[3][::-1], 0, 1) # reverse should not appear anywhere

print("\nTest Inexact - More Mismatch")
print("1 MATCH SHOULD OCCUR FOR EACH TEST")
print("-----------------------------------------")
print("insert")
test( dnas[1], 1, 1) # insertion
print("_______________")
print("delete")
test(dnas[1][0:3] + dnas[0][4:], 4, 2) # deletions
print("_______________")
print("substitute")
test(dnas[1][0:3] + 'T' + dnas[0][4:], 4, 2) # substitution

#input random protein
print("\nTest Inexact 2 - Combination of input")
#print(dnas[1])
print("_______________")
print("insert & del")
test( dnas[0][0:10] + 'C' + dnas[0][11:], 4, 2) # insert and delete
print("_______________")
print("insert & sub")
test(dnas[2][0:9] + 'T' + dnas[2][10:] + 'A',4,2) # insert and sub
print("_______________")
print("del & sub")
test(dnas[3][1:9] + 'T' + dnas[3][10:], 4, 2) # delete & substitution


#input random protein
print("\nTest Inexact 3 - All Combination of input")
#print(dnas[1])
print("_______________Test 4-1_______________")
test( dnas[4][0:10] + 'C' + dnas[4][11:15] + 'T' + dnas[4][16:] , 6, 2)
print("_______________Test 4-2_______________")
test( dnas[2][0:10] + 'C' + dnas[2][11:15] + 'A' + dnas[2][16:] , 6, 2)
print("_______________Test 4-3_______________")
test( dnas[1][0:10] + 'G' + dnas[1][11:15] + 'C' + dnas[1][16:] , 6, 2)


# Terminator Filter
print("\nTest Terminator Filter")
print("--------------------------------------")
print("Terminator in Middle of Alignment: No Alignment Output")
seq_t = "AKRSPQFTGKTMTQEERFEQR" # aligns to "AKRSPQFTGK&MTQEERFEQR"
dna_t = proteinToDNA(seq_t)
test(dna_t, 2, 2)
print("_______________")
print("Terminator at End of Alignment: Still Aligns")
seq_t = "REGVSAFLAKRSPQFTGKQM" # aligns to "REGVSAFLAKRSPQFTGK&M"
dna_t = proteinToDNA(seq_t)
dna_t = transcribe(dna_t)
test(dna_t, 4, 2)
print("_______________")
print("Terminator at Beginning of Alignment: Still Aligns")
seq_t = "KTMTQEERFEQRIAQETAIE" # aligns to "K&MTQEERFEQRIAQETAIE"
dna_t = proteinToDNA(seq_t)
dna_t = transcribe(dna_t)
test(dna_t, 4, 2)