import sys
from io import StringIO
import random 

# parse fasta, select random read, create fastq file 
input_text = sys.argv[1]
output_text = sys.argv[2]

# read in two strings from the test file 
def parse_fasta(fh): 
    reads = [] 
    while True: 
        first_line = fh.readline()
        if len(first_line) == 0:
            break 
        seq = first_line.rstrip()
        reads.append(seq)
    return reads

with open(input_text, 'r') as file:
    fasta = file.read()

result = parse_fasta(StringIO(fasta))

selected_reads = {}
# for the first five, reverse reads 
for i in range(10):
    read = ""
    idx = random.randint(0, 66281)
    for j in range(idx, idx+10):
        read += result[j]
    if i < 5:
        read = read[::-1]
    selected_reads[idx] = read

# generate fastq file using the above reads 
base = ""
for i in range(700): 
    base += "#"

with open(output_text, 'w') as f:
    for read_id, read in selected_reads.items():
        f.write("@r" + str(read_id))
        f.write("\n")
        f.write(read)
        f.write("\n")
        f.write("+")
        f.write("\n") 
        f.write(base)
        f.write("\n") 