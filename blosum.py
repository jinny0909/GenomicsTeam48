import numpy as np


def localAlignment(x, y, gap_default, gap_score, gap = True):


  BLOSUM62 = {
    "C":{"C":9, "S":-1, "T":-1, "P":-3, "A":0,  "G":-3, "N":-3, "D":-3, "E":-4, "Q":-3, "H":-3, "R":-3, "K":-3, "M":-1, "I":-1, "L":-1, "V":-1, "F":-2, "Y":-2, "W":-2},
    "S":{"C":-1,"S":4,  "T":1,  "P":-1, "A":1,  "G":0,  "N":1,  "D":0,  "E":0,  "Q":0,  "H":-1, "R":-1, "K":0,  "M":-1, "I":-2, "L":-2, "V":-2, "F":-2, "Y":-2, "W":-3},
    "T":{"C":-1,"S":1,  "T":4,  "P":1,  "A":-1, "G":1,  "N":0,  "D":1,  "E":0,  "Q":0,  "H":0,  "R":-1, "K":0,  "M":-1, "I":-2, "L":-2, "V":-2, "F":-2, "Y":-2, "W":-3},
    "P":{"C":-3,"S":-1, "T":1,  "P":7,  "A":-1, "G":-2, "N":-1, "D":-1, "E":-1, "Q":-1, "H":-2, "R":-2, "K":-1, "M":-2, "I":-3, "L":-3, "V":-2, "F":-4, "Y":-3, "W":-4},
    "A":{"C":0, "S":1,  "T":-1, "P":-1, "A":4,  "G":0,  "N":-1, "D":-2, "E":-1, "Q":-1, "H":-2, "R":-1, "K":-1, "M":-1, "I":-1, "L":-1, "V":-2, "F":-2, "Y":-2, "W":-3},
    "G":{"C":-3,"S":0,  "T":1,  "P":-2, "A":0,  "G":6,  "N":-2, "D":-1, "E":-2, "Q":-2, "H":-2, "R":-2, "K":-2, "M":-3, "I":-4, "L":-4, "V":0,  "F":-3, "Y":-3, "W":-2},
    "N":{"C":-3,"S":1,  "T":0,  "P":-2, "A":-2, "G":0,  "N":6,  "D":1,  "E":0,  "Q":0,  "H":-1, "R":0,  "K":0,  "M":-2, "I":-3, "L":-3, "V":-3, "F":-3, "Y":-2, "W":-4},
    "D":{"C":-3,"S":0,  "T":1,  "P":-1, "A":-2, "G":-1, "N":1,  "D":6,  "E":2,  "Q":0,  "H":-1, "R":-2, "K":-1, "M":-3, "I":-3, "L":-4, "V":-3, "F":-3, "Y":-3, "W":-4},
    "E":{"C":-4,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":2,  "E":5,  "Q":2,  "H":0,  "R":0,  "K":1,  "M":-2, "I":-3, "L":-3, "V":-3, "F":-3, "Y":-2, "W":-3},
    "Q":{"C":-3,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":0,  "E":2,  "Q":5,  "H":0,  "R":1,  "K":1,  "M":0,  "I":-3, "L":-2, "V":-2, "F":-3, "Y":-1, "W":-2},
    "H":{"C":-3,"S":-1, "T":0,  "P":-2, "A":-2, "G":-2, "N":1,  "D":1,  "E":0,  "Q":0,  "H":8,  "R":0,  "K":-1, "M":-2, "I":-3, "L":-3, "V":-2, "F":-1, "Y":2,  "W":-2},
    "R":{"C":-3,"S":-1, "T":-1, "P":-2, "A":-1, "G":-2, "N":0,  "D":-2, "E":0,  "Q":1,  "H":0,  "R":5,  "K":2,  "M":-1, "I":-3, "L":-2, "V":-3, "F":-3, "Y":-2, "W":-3},
    "K":{"C":-3,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":-1, "E":1,  "Q":1,  "H":-1, "R":2,  "K":5,  "M":-1, "I":-3, "L":-2, "V":-3, "F":-3, "Y":-2, "W":-3},
    "M":{"C":-1,"S":-1, "T":-1, "P":-2, "A":-1, "G":-3, "N":-2, "D":-3, "E":-2, "Q":0,  "H":-2, "R":-1, "K":-1, "M":5,  "I":1,  "L":2,  "V":-2, "F":0,  "Y":-1, "W":-1},
    "I":{"C":-1,"S":-2, "T":-2, "P":-3, "A":-1, "G":-4, "N":-3, "D":-3, "E":-3, "Q":-3, "H":-3, "R":-3, "K":-3, "M":1,  "I":4,  "L":2,  "V":1,  "F":0,  "Y":-1, "W":-3},
    "L":{"C":-1,"S":-2, "T":-2, "P":-3, "A":-1, "G":-4, "N":-3, "D":-4, "E":-3, "Q":-2, "H":-3, "R":-2, "K":-2, "M":2,  "I":2,  "L":4,  "V":3,  "F":0,  "Y":-1, "W":-2},
    "V":{"C":-1,"S":-2, "T":-2, "P":-2, "A":0,  "G":-3, "N":-3, "D":-3, "E":-2, "Q":-2, "H":-3, "R":-3, "K":-2, "M":1,  "I":3,  "L":1,  "V":4,  "F":-1, "Y":-1, "W":-3},
    "F":{"C":-2,"S":-2, "T":-2, "P":-4, "A":-2, "G":-3, "N":-3, "D":-3, "E":-3, "Q":-3, "H":-1, "R":-3, "K":-3, "M":0,  "I":0,  "L":0,  "V":-1, "F":6,  "Y":3,  "W":1},
    "Y":{"C":-2,"S":-2, "T":-2, "P":-3, "A":-2, "G":-3, "N":-2, "D":-3, "E":-2, "Q":-1, "H":2,  "R":-2, "K":-2, "M":-1, "I":-1, "L":-1, "V":-1, "F":3,  "Y":7,  "W":2},
    "W":{"C":-2,"S":-3, "T":-3, "P":-4, "A":-3, "G":-2, "N":-4, "D":-4, "E":-3, "Q":-2, "H":-2, "R":-3, "K":-3, "M":-1, "I":-3, "L":-2, "V":-3, "F":1,  "Y":2,  "W":11}}
  
  perfect_score = 0 
  for c in y: 
    perfect_score += BLOSUM62[c][c]
  
  V = np.zeros((len(x)+1, len(y)+1), dtype=int)
  for i in range(1, len(x)+1):
    for j in range(1, len(y)+1):
      V[i, j] = max(V[i-1, j-1] + BLOSUM62[x[i-1]][y[j-1]], # diagonal
      V[i-1, j] + gap_score, # vertical
      V[i, j-1] + gap_score, # horizontal
      0) # empty
      argmax = np.where(V == V.max())

  print("PERFECT alignment score should be: \n")
  print(perfect_score)

  print("Actual alignment score: \n")
  print(int(V[argmax]))
  return int(V[argmax])




# main : 
# call the above function as follows: 

# localAlignment(ref, pattern, gap_default, gap_score, gap = True)