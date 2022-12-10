# Computational Genomics Final Project
# Team 48

import time
from ctypes import alignment
from bwa import *
from proteinToDNA import *
from proteinIndex import *
from openRF import *
from blosum import *
import random
from parse_data import *

""" Full implementation of tool in a single read format"""
# mismatch : maximum mismatch threshold
# check : 0 for exact matching with FM Index validation
def align_tool(fasta_file, read, mismatch, check):
    # parse database
    reference, proteinID, proteins = parse_data(fasta_file)
    SA, C, O, O_rev, n, alphabet, proteinIDs = precalculation(reference + '$')
    protein_index = FMIndex(reference + '$', 10, 5) # optional for exact matching

    # get reading frames
    tests = getRF(read)
    passed_test, passed_protein = [], []

    # loop through possible frames
    for test1 in tests:
        # call bw alignment
        output = bwa(test1, SA, C, O, O_rev, n, alphabet, proteinIDs, mismatch)
        # check found matches
        for o in output:
            passed_test.append(test1)
            passed_protein.append(proteins[proteinID[o]])
            # output
            print("Protein match: %s" % proteinID[o])

            # check alignment
            perfect_score, alignment_score = localAlignment(proteins[proteinID[o]], test1, -1, -1, True)

            print("Perfect Alignment Score:")
            print(perfect_score)
            print('Alignment Score')
            print(alignment_score)
            print("\n")

    if check == 0:
        for p in passed_test:
            # use FM Index to check
            output = protein_index.query(p)
            print(output) # this is index in reference frame that matches
            if len(output) == 1:
                print("Validated: Alignment Found")
