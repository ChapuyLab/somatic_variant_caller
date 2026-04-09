from pathlib import Path
import tempfile
import warnings

"""
Setting up variables for work, log and scratch_dir
"""
if "workdir" in config:
    wrkdir = Path(config["workdir"])
else:
    raise ValueError("Please define a workdir")

if "logdir" in config:
    logdir = Path(config["logdir"])
else:
    raise ValueError("Please define a logdir")

# if "scratch_dir" not in config:
scratch_dir = tempfile.gettempdir()
# else:
#     scratch_dir = config[
#         "scratch_dir"
#     ]  # snakemake doesnt accept path in resources has to be a string
#     if not os.path.exists(scratch_dir):
#         print(get_data_time(), "Scratch directory does not exist")
#         os.makedirs(scratch_dir)

"""
Setting up variables related to sample information
"""

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

samples = list()  # kept for back compatibility

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


"""
Setup variables for database used by mutect2
"""
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

"""
The target file only to declared for WES and Panel data
Scatter count that should be used
Only use for WES or Panel sequencing otherwise
please download the scatter list from
https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37/scattered_wgs_intervals/scatter-80
which divides the genome into 80 intervals
"""

if "target_file" in config:
    print("Using Targeted/WES Sequencing mode")
    target_file = config["target_file"]
    scatter_count = int(config["scatter_count"]) if "scatter_count" in config else 10
    interval_folder = False
    # beware please make the scatter appropriate to the number of intervals you have in bed file#
else:
    print("Using WGS mode")
    if "interval_folder" in config:
        scatter_count = 80
        interval_folder = Path(config["interval_folder"])
    else:
        print("Interval folder not found in config")
        scatter_count = 22
    target_file = False

"""
Set up variables related to annovar
"""
run_annovar = config["run_annovar"] if "run_annovar" in config else False

if run_annovar:
    path_to_annovar = Path(config["path_to_annovar"])
    annovar_db = Path(config["annovar_db"])

    if "protocol" in config and "operation" in config and "geneDatabase" in config:
        annovar_protocol = config["protocol"]
        annovar_operation = config["operation"]
        geneDatabase = config["geneDatabase"]
    else:
        if "protocol" not in config or "operation" in config:
            print(
                "Both Protocol and operation need to be defined for annovar, falling back to default"
            )
        annovar_protocol = (
            "wgEncodeGencodeBasicV19,cytoBand,avsnp150,dbnsfp42c,clinvar_20240917"
        )
        annovar_operation = "g,r,f,f,f"
        geneDatabase = "wgEncodeGencodeBasicV19"


"""
Set up variables for Sample swap check betweem tumor and normal samples
"""
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

"""
Setup variables for detin
"""
if "pkl_detin" not in config:
    raise ValueError("Please provide pickle file")
else:
    pickle_high_af = config["pkl_detin"]

"""
Set the variable for blat filter
"""
if "genome_build" not in config:
    raise ValueError("Please provide genome build")
else:
    genome_build = config["genome_build"]

if "insertSize" in config:
    insertSize = config["insertSize"]
else:
    insertSize = None

stepper = config["stepper"] if "stepper" in config else None

if "fake_normal" in config and config["fake_normal"]:
    print("Running with a fake normal")
    fake_normal = config["fake_normal"]
else:
    fake_normal = False
