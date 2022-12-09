# Computational Genomics Final Project
# Team 48
"""This file is for testing purposes.  It can be run to replicate the results in the writeup, and it uses all
 the methods created for this project.  The three tests are for exact protein alignment, inexact protein alignment
 and checking that terminators do not match."""
import time
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
protein_index = FMIndex(reference + '$', 10, 5)

# output = bwa("RRYPHREPVPHRVQHLQRLRHHLRPYPVPRQQPLC", SA, C, O, O_rev, n, alphabet, proteinIDs, 4)

def test(DNA, mismatch, check):
    tests = getRF(DNA)
    passed_test, passed_protein = [], []

    for test1 in tests:
        start = time.time()
        output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs, mismatch)
        end = time.time()
        print('BWA elapsed time (millisec):')
        print((end-start)*1000)
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
    # else:
        # for test1 in tests:
        #     output = protein_index.query(test1)
        #     assert(output == (-1, -1))

name_list = [' 3-isopropylmalate dehydrogenase ', ' 2-iminoacetate synthase ', ' 1,4-dihydroxy-2-naphthoyl-CoA synthase ', 
' 2,5-diketo-D-gluconic acid reductase B ', ' 2-hydroxy-3-oxopropionate reductase ' ]

#name_list = [' 2,5-diketo-D-gluconic acid reductase A ',' 2-keto-3-deoxygluconate permease ', ' 3-deoxy-D-manno-octulosonate 8-phosphate phosphatase KdsC ',
#' 3-mercaptopyruvate sulfurtransferase ', ' 3-oxoacyl-[acyl-carrier-protein] synthase 1 ' ]
dnas = []
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

# Check that mismatch gives no Exact Matches
# print("NO MATCHES SHOULD OCCUR")
print("-----------------------------------------")
test(dnas[0][0:10] + 'C' + dnas[0][10:], 0, 1) # insertion
test(dnas[1][0:9] + dnas[1][10:], 0, 1) # deletions
test(dnas[2][0:9] + 'T' + dnas[2][10:], 0, 1) # substitution
test(dnas[3][::-1], 0, 1) # reverse should not appear anywhere

print("Test Inexact 2 - More Mismatch")
#print(dnas[1])
print("_______________")
print("insert")
test( dnas[1], 1, 1) # insertion
print("_______________")
print("delete")
test(dnas[1][0:3] + dnas[0][4:], 4, 1) # deletions
print("_______________")
print("substitute")
test(dnas[1][0:3] + 'T' + dnas[0][4:], 4, 1) # substitution

#input random protein
print("Test Inexact 3 - Combination of input")
#print(dnas[1])
print("_______________")
print("insert & del")
test( dnas[0][0:10] + 'C' + dnas[0][11:], 4, 1) # insert and delete
print("_______________")
print("insert & sub")
test(dnas[2][0:9] + 'T' + dnas[2][10:] + 'A',4,1) # insert and sub
print("_______________")
print("del & sub")
test(dnas[3][1:9] + 'T' + dnas[3][10:], 4, 1) # delete & substitution


#input random protein
print("Test Inexact 4 - All Combination of input")
#print(dnas[1])
print("_______________Test 4-1_______________")
test( dnas[4][0:10] + 'C' + dnas[4][11:15] + 'T' + dnas[4][16:] , 6, 1) 
print("_______________Test 4-2_______________")
test( dnas[2][0:10] + 'C' + dnas[2][11:15] + 'A' + dnas[2][16:] , 6, 1) 
print("_______________Test 4-3_______________")
test( dnas[1][0:10] + 'G' + dnas[1][11:15] + 'C' + dnas[1][16:] , 6, 1) 
