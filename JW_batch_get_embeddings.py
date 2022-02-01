# def add_hash(x,y):
#   return x+"_"+hashlib.sha1(y.encode()).hexdigest()[:5]

# query_sequence = 'PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK'
#
# jobname = 'test'
# # remove whitespaces
# basejobname = "".join(jobname.split())
# basejobname = re.sub(r'\W+', '', basejobname)
# jobname = add_hash(basejobname, query_sequence)
# while os.path.isfile(f"{jobname}.csv"):
#   jobname = add_hash(basejobname, ''.join(random.sample(query_sequence,len(query_sequence))))
#
# with open(f"{jobname}.csv", "w") as text_file:
#   text_file.write(f"id,sequence\n{jobname},{query_sequence}")
jobname = 'test_3a301'
queries_path = f"{jobname}.csv"

use_amber = False
use_templates = False
save_to_google_drive = False

msa_mode = "MMseqs2 (UniRef+Environmental)" #@param ["MMseqs2 (UniRef+Environmental)", "MMseqs2 (UniRef only)","single_sequence","custom"]
model_type = "auto"
pair_mode = "unpaired+paired"
num_recycles = 1

# # decide which a3m to use
# if msa_mode.startswith("MMseqs2"):
#   a3m_file = f"{jobname}.a3m"
# elif msa_mode == "custom":
#     pass
# else:
#   a3m_file = f"{jobname}.single_sequence.a3m"
#   with open(a3m_file, "w") as text_file:
#     text_file.write(">1\n%s" % query_sequence)


from pathlib import Path

from colabfold.batch import get_queries, run_get_embeddings, set_model_type
from colabfold.download import download_alphafold_params
from colabfold.utils import setup_logging

result_dir="."
setup_logging(Path(".").joinpath("log.txt"))
queries, is_complex = get_queries(queries_path)
model_type = set_model_type(is_complex, model_type)
download_alphafold_params(model_type, Path("."))

predictions_list = run_get_embeddings(
    queries=queries,
    result_dir=result_dir,
    use_templates=use_templates,
    use_amber=use_amber,
    msa_mode=msa_mode,
    model_type=model_type,
    num_models= 1,#5,
    num_recycles=num_recycles,
    model_order=[2],#[1, 2, 3, 4, 5],
    is_complex=is_complex,
    data_dir=Path("."),
    keep_existing_results=False,
    recompile_padding=1.0,
    rank_by="auto",
    pair_mode=pair_mode,
    stop_at_score=float(100),
)
first_amino_acid_embedding = predictions_list[0]['representations']['msa_first_row'][0]
breakpoint_var = True