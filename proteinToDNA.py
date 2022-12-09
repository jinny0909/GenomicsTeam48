
def reverseComplement(seq): 
  complement = {"A": "U", "U": "A", "T": "A", "C": "G", "G": "C"}
  output = ''
  for c in seq: 
    if c in 'ACGT':
      output += complement[c]

  return output


def proteinToDNA(seq): 
    complement = {"U": "A", "A": "T", "G" : "C", "C":"G"}
    lookup = {"UUU": "F",    
            "CUU": "L", "AUU": "I", "GUU": "V",
            "UUC": "F",    
            "CUC": "L", "AUC": "I", 
            "GUC": "V", "UUA": "L", "CUA": "L", "AUA": "I", 
            "GUA": "V","UUG": "L",    
            "CUG": "L", "AUG": "M", "GUG": "V","UCU": "S",    
            "CCU": "P", "ACU": "T", "GCU": "A","UCC": "S",    
            "CCC": "P", "ACC": "T", "GCC": "A","UCA": "S",    
            "CCA": "P", "ACA": "T", "GCA": "A","UCG": "S",    
            "CCG": "P", "ACG": "T", "GCG": "A","UAU": "Y",    
            "CAU": "H", "AAU": "N", "GAU": "D","UAC": "Y",    
            "CAC": "H", "AAC": "N", "GAC": "D","UAA": "Stop", 
            "CAA": "Q", "AAA": "K", "GAA": "E","UAG": "Stop", 
            "CAG": "Q", "AAG": "K", "GAG": "E","UGU": "C",    
            "CGU": "R", "AGU": "S", "GGU": "G","UGC": "C",    
            "CGC": "R", "AGC": "S", "GGC": "G","UGA": "Stop", 
            "CGA": "R", "AGA": "R", "GGA": "G","UGG": "W",    
            "CGG": "R", "AGG": "R", "GGG": "G"}


    rev_lookup = {'F': 'UUC', 'L': 'CUG', 'I': 'AUA', 'V': 'GUG',
     'M': 'AUG', 'S': 'AGC', 'P': 'CCG', 'T': 'ACG', 'A': 'GCG', 
     'Y': 'UAC', 'H': 'CAC', 'N': 'AAC', 'D': 'GAC', 'Stop': 'UGA', 
     'Q': 'CAG', 'K': 'AAG', 'E': 'GAG', 'C': 'UGC', 'R': 'AGG', 
     'G': 'GGG', 'W': 'UGG'}
   

    text = ""
    for c in seq: # for each char in seq
        text += rev_lookup[c]
    
    output = ""
    for c in text: 
        output += complement[c]

    return output

# main: call this function with the given format 
# proteinToDNA("FLIMVSPTAYHQNKDECWRG")