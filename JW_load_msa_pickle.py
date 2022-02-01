import os
import pickle

with open('saved_msa/1SW6B.pickle', 'rb') as f:
    msa = pickle.load(f)

breakpoint_var = True

with open('saved_msa/examined_msa.txt', 'w') as f2:
    f2.write(msa['unpaired_msa'][0])