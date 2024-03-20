from pathlib import Path

normal = config["normal_sample"]
tumour = config["tumour_sample"]

wrkdir = Path(config["workdir"])


logdir = Path(config["logdir"])

metadata = config["metadata"]

samples = [normal, tumour]
common_biallelic = config["common_biallelic"]
germline = config["germline"]
genome = config["genome"]

target_file = config["target_file"] if "target_file" in config else False

scatter_count = 10  # beware please make the scatter appropriate to the number of intervals you have in bed file#

run_annovar = config["run_annovar"] if "run_annovar" in config else False

path_to_annovar = Path(config["path_to_annovar"])
annovar_db = Path(config["annovar_db"])


rg_normal = config["rg_normal"]
rg_tumour = config["rg_tumour"]
