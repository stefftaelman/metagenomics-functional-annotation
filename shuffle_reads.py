from Bio import Seq, SeqIO, SeqRecord
import numpy as np
import os
import pandas as pd
import random
from tqdm import tqdm

# BEST TO RUN THIS ON A VM

READ_DIR = "output-mount/fasta-files/"
OUT_DIR = "output-mount/shuffled-reads/"
N = 10

# randomly pick 10 samples whose reads to shuffle
list_of_reads = []
for root, dirs, files in os.walk(READ_DIR):
	for file in files:
         if ".DS_Store" not in file:
            list_of_reads.append(os.path.join(root,file))
to_shuffle_reads = np.random.choice(list_of_reads, N)
print(f"{len(to_shuffle_reads)} reads were chosen at random to shuffle the contents of their reads.")


# loop over the to-shuffle-reads and randomize their contents (preserving nucleotide fractions) 
for file in tqdm(to_shuffle_reads):
    records = []
    fname = OUT_DIR + file.split('/')[-1].split('m')[0] + "shuffled.fasta"
    fasta_sequences = SeqIO.parse(open(file),'fasta')
    for seq in fasta_sequences:
        seq_header = seq.description + " | SHUFFLED"
        read = seq.seq
        shuffled_read = ''.join(list(random.sample(list(read), len(read))))
        seq_obj = Seq.Seq(shuffled_read)
        records.append(SeqRecord.SeqRecord(seq_obj, seq_header))
    SeqIO.write(records, fname, "fasta")
    print(f"Wrote out shuffled reads to {fname}!")