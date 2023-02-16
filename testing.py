from Bio import SeqIO
import os

print(os.getcwd())
fasta_sequences = SeqIO.parse(open("coursea_course/res/datasets/upstream250.txt"), 
'fasta')

print(list(fasta_sequences))

import sys
print(sys.executable)
#%%
