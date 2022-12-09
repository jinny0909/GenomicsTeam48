# This file contains functions related to creating reading frames 
def reverseComplement(seq): 
  complement = {"A": "U", "U": "A", "T": "A", "C": "G", "G": "C"}
  output = ''
  for c in seq: 
    if c in 'ACGT':
      output += complement[c]

  return output

def transcribe(seq):
  output = ""
  for c in seq:
    if c in "ACG":
      output += c
    elif c == "T":
      output += "U"

def translate(seq):
  lookup = {
  "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
  "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
  "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
  "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
  "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
  "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
  "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
  "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
  "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
  "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
  "UAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
  "UAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
  "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
  "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
  "UGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
  "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
  codon, protein = "", ""
  for c in seq:
    if c in "ACGU":
      codon += c
    if len(codon) == 3:
      protein += lookup[codon]
      codon = ""

  return protein

def longestORF(seq):
  longest = 0
  current = 0
  for c in seq:
    if c == "*":
      if current > longest:
        longest = current
      current = 0
    else:
      current+=1
  return longest

def bestFrame(seqs):
  best = []
  longestORF = 0
  for s in seqs:
    ORF = longestORF(s)
    if ORF > longestORF:
      best = [s]
    elif ORF == longestORF:
      best.append([s])
  return best