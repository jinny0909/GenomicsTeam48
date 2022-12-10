# Computational Genomics
# Team 48
"""This function parses a fasta file with the protein sequences and returns the concatenated reference string,
the list of protein names, and a dataset with the protein sequences themselves"""
def parse_data(dataset):
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

    return reference, protein_id, proteins