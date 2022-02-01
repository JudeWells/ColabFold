import pickle
from colabfold.batch import get_msa_and_templates, msa_to_str, generate_input_feature, predict_structure_get_embeddings, load_models_and_params
from pathlib import Path
import pandas as pd
import os
# queries is a list which is passed into the outer function above this which is run_get_embeddings()
# [('test_7dfa6', 'PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK', None)]

def pickler(jobname, unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features):
    msa_dict = {
        'unpaired_msa': unpaired_msa,
        'paired_msa': paired_msa,
        'query_seqs_unique': query_seqs_unique,
        'query_seqs_cardinality': query_seqs_cardinality,
        'template_features': template_features,

    }
    directory_out = 'saved_msa/'
    filename = directory_out + jobname + '.pickle'
    with open(filename, 'wb') as handle:
        pickle.dump(msa_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Set up configuration and settings
msa_mode = "MMseqs2 (UniRef+Environmental)" #  "MMseqs2 (UniRef only)"
result_dir = Path('')
use_templates = False
pair_mode = "unpaired+paired"
host_url = 'https://a3m.mmseqs.com'

# load the dataframe which contains the sequences
filepath = '/Users/judewells/Documents/dataScienceProgramming/cath-funsite-predictor/experiments/PPI_training_dataset_with_sequences.csv'
df = pd.read_csv(filepath)
df['domain_id'] = (df.PDBID + df.CHAIN).str.upper()

# generate the MSAs for each chain
for i, row in df.iterrows():
    query_sequence = row.sequence
    if query_sequence is None:
        continue
    jobname = row.domain_id

    msa_results = get_msa_and_templates(jobname, query_sequence, result_dir,
                    msa_mode, use_templates, pair_mode, host_url)

    # split the msa_results tuple into its constitutent parts
    unpaired_msa,paired_msa,query_seqs_unique,query_seqs_cardinality,template_features  = msa_results
    pickler(jobname, unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features)
    if i % 50 == 0:
        print(f'completed {i}')










###### The rest of the code is not essential for the allignent
###### but shows which aspects of it are required

msa = msa_to_str(
    unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality
)

is_complex = False
use_amber = False
model_type = 'AlphaFold2-ptm'
crop_len = 0
rank_by = 'plddt'
stop_at_score = 100
prediction_callback = None
num_models = 1
num_recycles = 1
model_order = [2]
model_extension = '_ptm'
data_dir = Path('')
recompile_all_models = False





query_sequence_len_array = [
    len(query_seqs_unique[i])
    for i, cardinality in enumerate(query_seqs_cardinality)
    for _ in range(0, cardinality)
]

model_runner_and_params = load_models_and_params(
    num_models,
    use_templates,
    num_recycles,
    model_order,
    model_extension,
    data_dir,
    recompile_all_models,
    stop_at_score=stop_at_score,
    rank_by=rank_by,
)


input_features = generate_input_feature(
    query_seqs_unique,
    query_seqs_cardinality,
    unpaired_msa,
    paired_msa,
    template_features,
    is_complex,
    model_type,
)

predictions_list = predict_structure_get_embeddings(
    jobname,
    result_dir,
    input_features,
    is_complex,
    use_templates,
    sequences_lengths=query_sequence_len_array,
    crop_len=crop_len,
    model_type=model_type,
    model_runner_and_params=model_runner_and_params,
    do_relax=use_amber,
    rank_by=rank_by,
    stop_at_score=stop_at_score,
    prediction_callback=prediction_callback,
)

