import pickle
from colabfold.batch import get_msa_and_templates, msa_to_str, generate_input_feature, predict_structure_get_embeddings, load_models_and_params
from pathlib import Path
import pandas as pd
import os
import time
# queries is a list which is passed into the outer function above this which is run_get_embeddings()
# [('test_7dfa6', 'PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK', None)]

def pickler(jobname, unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features, directory_out = 'saved_msa/'):
    chain_ids = jobname.split('-|-')
    msa_dict = {
        'unpaired_msa': unpaired_msa,
        'paired_msa': paired_msa,
        'query_seqs_unique': query_seqs_unique,
        'query_seqs_cardinality': query_seqs_cardinality,
        'template_features': template_features,

    }
    filename = directory_out + jobname + '.pickle'
    with open(filename, 'wb') as handle:
        pickle.dump(msa_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Set up configuration and settings
msa_mode = "MMseqs2 (UniRef+Environmental)" #  "MMseqs2 (UniRef only)"
result_dir = Path('')
use_templates = False
pair_mode = "unpaired+paired"
host_url = 'https://a3m.mmseqs.com'
directory_out = 'saved_msa/'

# load the dataframe which contains the sequences
filepath = '/Users/judewells/Documents/dataScienceProgramming/cath-funsite-predictor/experiments/PPI_training_dataset_with_sequences.csv'
df = pd.read_csv(filepath)
df['domain_id'] = (df.PDBID + df.CHAIN).str.upper()

# generate the MSAs for each chain
# for i, row in df.iterrows():
#     query_sequence = row.sequence
#     if query_sequence is None:
#         continue
#     jobname = row.domain_id
#     filepath = directory_out + jobname + '.pickle'
#     if os.path.exists(filepath):
#         continue

jobname = '-|-'.join(df.loc[70:74,'domain_id'].values)
batch_of_sequences = df.loc[70:74,'sequence'].values

msa_results = get_msa_and_templates(jobname, batch_of_sequences, result_dir,
                msa_mode, use_templates, pair_mode, host_url)

# split the msa_results tuple into its constitutent parts
unpaired_msa,paired_msa,query_seqs_unique,query_seqs_cardinality,template_features  = msa_results
pickler(jobname, unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality,
        template_features, directory_out = 'saved_msa/')




