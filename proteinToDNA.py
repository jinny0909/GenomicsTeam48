def proteinToDNA(seq):
    rev_lookup = {'F': 'UUC', 'L': 'CUG', 'I': 'AUA', 'V': 'GUG',
     'M': 'AUG', 'S': 'AGC', 'P': 'CCG', 'T': 'ACG', 'A': 'GCG', 
     'Y': 'UAC', 'H': 'CAC', 'N': 'AAC', 'D': 'GAC', '*': 'UGA',
     'Q': 'CAG', 'K': 'AAG', 'E': 'GAG', 'C': 'UGC', 'R': 'AGG', 
     'G': 'GGG', 'W': 'UGG'}
   

    text = ""
    for c in seq: # for each char in seq
        text += rev_lookup[c]

    return text

# main: call this function with the given format 
# proteinToDNA("FLIMVSPTAYHQNKDECWRG")