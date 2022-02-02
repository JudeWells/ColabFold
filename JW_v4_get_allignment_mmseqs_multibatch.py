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
    for idx, ch_id in enumerate(chain_ids):
        msa_dict = {
            'unpaired_msa': unpaired_msa[idx],
            'paired_msa': paired_msa[idx],
            'query_seqs_unique': query_seqs_unique[idx],
            'query_seqs_cardinality': query_seqs_cardinality[idx],
            'template_features': template_features[idx],

        }
        filename = directory_out + ch_id + '.pickle'
        with open(filename, 'wb') as handle:
            pickle.dump(msa_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Set up configuration and settings
msa_mode = "MMseqs2 (UniRef+Environmental)" #  "MMseqs2 (UniRef only)"
result_dir_name = 'result_dir'
if not os.path.exists(result_dir_name):
    os.mkdir(result_dir_name)
result_dir = Path(result_dir_name)
use_templates = False
pair_mode = "unpaired+paired"
host_url = 'https://a3m.mmseqs.com'
directory_out = 'saved_msa/'
if not os.path.exists(directory_out):
    os.mkdir(directory_out)
# load the dataframe which contains the sequences
filepath = '../data/training_data.csv'
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
n_seqs_per_batch = 5

for i in range(1000,df.last_valid_index(), n_seqs_per_batch):
    last_index = i+n_seqs_per_batch-1
    if last_index > df.last_valid_index():
        last_index = df.last_valid_index()
    one_batch = df.loc[i:last_index, :]
    chain_ids = one_batch.domain_id.values
    batch_completed = False
    for chid in chain_ids:
        if chid + '.pickle' in os.listdir(directory_out):
            batch_completed = True
    if not batch_completed:
        print(f'\n\n{i} sequences complete\n\n')
        jobname = '-|-'.join(chain_ids)
        batch_of_sequences = one_batch.sequence.values
        try:
            msa_results = get_msa_and_templates(jobname, batch_of_sequences, result_dir,
                            msa_mode, use_templates, pair_mode, host_url)

            # split the msa_results tuple into its constitutent parts
            unpaired_msa,paired_msa,query_seqs_unique,query_seqs_cardinality,template_features = msa_results
            pickler(jobname, unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality,
                    template_features, directory_out='saved_msa/')

        except:
            pass


