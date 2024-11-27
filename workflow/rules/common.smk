from pathlib import Path
import tempfile


samples = list()
if "normal_sample" in config:
    normal = config["normal_sample"]
    samples.append(normal)
else:
    normal = None
    print("Using tumour only mode")

if "tumour_sample" in config:
    tumour = config["tumour_sample"]
    samples.append(tumour)
else:
    raise ValueError("Please define a tumour sample")

if "workdir" in config:
    wrkdir = Path(config["workdir"])
else:
    raise ValueError("Please define a workdir")

if "logdir" in config:
    logdir = Path(config["logdir"])
else:
    raise ValueError("Please define a logdir")

if "scratch_dir" not in config:
    scratch_dir = tempfile.gettempdir()
else:
    scratch_dir = config[
        "scratch_dir"
    ]  # snakemake doesnt accept path in resources has to be a string
    if not os.path.exists(scratch_dir):
        print(get_data_time(), "Scratch directory does not exist")
        os.makedirs(scratch_dir)

if "metadata" not in config:
    raise ValueError("Please define a metadata file, or bam file paths")
else:
    metadata = config["metadata"]

if "normal_bam_file" in config:
    normal_bam_file = config["normal_bam_file"]
else:
    normal_bam_file = False

if "tumour_bam_file" in config:
    tumour_bam_file = config["tumour_bam_file"]
else:
    tumour_bam_file = False


if "pid" in config:
    pid = config["pid"]
else:
    raise ValueError("Please define a pid")


if "common_biallelic" not in config:
    raise ValueError("Please define a common_biallelic file")
else:
    common_biallelic = config["common_biallelic"]

if "germline" not in config:
    raise ValueError("Please define a gnomad file")
else:
    germline = config["germline"]

if "genome" not in config:
    raise ValueError("Please provide a genome file in config")
else:
    genome = config["genome"]

if "pon" not in config:
    print("Pipeline will be run without a panel of normals")
    pon = False
else:
    pon = config["pon"]


target_file = config["target_file"] if "target_file" in config else False


scatter_count = (
    int(config["scatter_count"]) if "scatter_count" in config else 22
)  # beware please make the scatter appropriate to the number of intervals you have in bed file#


run_annovar = config["run_annovar"] if "run_annovar" in config else False

path_to_annovar = Path(config["path_to_annovar"])
annovar_db = Path(config["annovar_db"])

if "rg_normal" not in config:
    rg_normal = None
else:
    rg_normal = config["rg_normal"]

if "rg_tumour" not in config:
    rg_tumour = None
else:
    rg_tumour = config["rg_tumour"]


if normal is not None:
    output_prefix = str(pid) + "_" + tumour + "_" + normal
else:
    output_prefix = str(pid) + "_" + tumour

if "LOD" not in config:
    raise ValueError("Please define a LOD file")
else:
    LOD = config["LOD"]

if "EXPECT_ALL_GROUPS_TO_MATCH" not in config:
    raise ValueError("Please define a EXPECT_ALL_GROUPS_TO_MATCH")
else:
    expect_all_groups_to_match = config["EXPECT_ALL_GROUPS_TO_MATCH"]

if "haplotypemap" not in config:
    raise ValueError("Please define a haplotypemap file")
else:
    haplotypemap = config["haplotypemap"]

if "pkl_detin" not in config:
    raise ValueError("Please provide pickle file")
else:
    pickle_high_af = config["pkl_detin"]

if "genome_build" not in config:
    raise ValueError("Please provide genome build")
else:
    genome_build = config["genome_build"]

if "insertSize" in config:
    insertSize = config["insertSize"]
else:
    insertSize = None
# blat_binary = "/dh-projects/ag-ishaque/analysis/sahays/pipelines/mutect2_pipeline/executables/blat"
# faToTwoBit_binary = "/dh-projects/ag-ishaque/analysis/sahays/pipelines/mutect2_pipeline/executables/faToTwoBit"
stepper = "all"
