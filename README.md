# ChapuyLab Somatic Variant Caller

This repository contains a Snakemake-based pipeline for running somatic variant calling using GATK Mutect2. It is designed to process individual samples (both tumor-normal paired and tumor-only) efficiently on HPC clusters.

The setup and job orchestration are managed via a helper Jupyter Notebook, which generates sample-specific configurations and execution scripts.

---

## Features

* **Sample Configurations:** Supports matched Tumor-Normal pairs, Tumor-Only runs, and "Fake Normal" fallbacks.
* **Cluster Ready:** Pre-configured for SLURM execution, but easily adaptable to other scheduling engines via Snakemake profiles.
* **Annotation & Filtering:** Integrates ANNOVAR for variant annotation and deTiN for estimating tumor-in-normal contamination.
* **Containerized:** Leverages Conda/Mamba and Singularity for reproducible environments.
* **Tested Build:** The pipeline has been extensively tested on the `hg19` genome build.

---

## Installation

### 1. Clone the Repository
First, clone the pipeline to your local environment or HPC cluster and navigate into the directory:
```bash
git clone [https://github.com/ChapuyLab/somatic_variant_caller.git](https://github.com/ChapuyLab/somatic_variant_caller.git)
cd somatic_variant_caller
```

### 2. Base Environment Setup
The pipeline relies on a Conda/Mamba environment to run correctly. The `base_env` configuration files are present at the root path of the cloned repository. Create and activate the environment before proceeding:
```bash
# Example assuming the file is named environment.yml or base_env.yml
mamba env create -f base_env.yml
mamba activate base_env
```

---

## Prerequisites and Required Data

Before running the setup notebook, you must download and configure several reference databases and tools.

### 1. Reference Files (GATK Best Practices)
Download the following from the Google Cloud GATK Best Practices bucket:
* **Germline resource:** `af-only-gnomad.raw.sites.vcf`
* **Common biallelic variants:** `small_exac_common_3.vcf`
* **Panel of Normals:** `pon.vcf`
* **Reference Genome:** `hs37d5_PhiX.fa` (or the standard Broad `hg19` build)
* **Haplotype map:** `hg19_nochr.map` (Available via [fingerprint_maps](https://github.com/naumanjaved/fingerprint_maps/tree/master/map_files))

### 2. Annotation & Filtering Resources
* **deTiN Data:** Requires high AF ExAC pickle file (`exac.pickle_high_af`).
* **Target Baits:** Please provide a target or baits file (`Agilent_SureSelectv2_baits.bed`). Baits is recommended as it helps with downstream cnvcalling

* **ANNOVAR:** You must download the ANNOVAR software and its human database (`humandb`). *See the detailed setup guide below.*

---

## ANNOVAR Setup and Database Installation Guide

This section provides step-by-step instructions for installing ANNOVAR and downloading the specific databases required for the variant annotation protocol: `wgEncodeGencodeBasicV19`, `cytoBand`, `avsnp150`, `dbnsfp42c`, and `clinvar_20240917`.


### Prerequisites
* A Linux/Unix or macOS environment.
* **Perl** installed on your system (standard on most Unix systems).
* Internet connection to download the software and databases.

### Step 1: Download and Extract ANNOVAR

1.  **Register and Download**: ANNOVAR is free for non-commercial use, but it requires registration. Go to the [ANNOVAR website](https://annovar.openbioinformatics.org/) and fill out the registration form. You will receive an email with a link to download the `annovar.latest.tar.gz` package.
2.  **Extract the Package**: Once downloaded, move the file to your desired installation directory and run the following command in your terminal:
    ```bash
    tar -zxvf annovar.latest.tar.gz
    ```
3.  **Navigate to the Directory**:
    ```bash
    cd annovar/
    ```
    *You should see several Perl scripts here, including `annotate_variation.pl` and `table_annovar.pl`, as well as a folder named `humandb/`.*

### Step 2: Download the Required Databases

You will use the `annotate_variation.pl` script to pull down your specific databases into the `humandb/` directory. Run the following commands one by one.

> **Storage Warning:** Some of these databases (like dbNSFP) are very large. Ensure you have sufficient disk space (at least 20-30 GB free) and a stable internet connection.

```bash
# 1. Gene-based annotation (GENCODE V19)
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar wgEncodeGencodeBasicV19 humandb/

# 2. Region-based annotation (Cytogenetic bands)
perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/

# 3. Filter-based annotation (dbSNP 150)
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/

# 4. Filter-based annotation (dbNSFP v4.2c - Functional predictions)
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42c humandb/

# 5. Filter-based annotation (ClinVar - Sept 17, 2024)
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20240917 humandb/
```
### 3. Configure the Snakemake Profile
Before running the setup notebook or any jobs, you **must** edit the Snakemake profile configuration file located at `profile/config.yaml` to match your cluster's architecture.

Open the file and modify the following critical parameters:

* **`conda-prefix` and `singularity-prefix`**: These **must** be changed to absolute paths pointing to directories where you want to store your environments and images (e.g., `/home/user/somatic_variant_caller/conda_envs`). If you leave them as relative paths, Snakemake will redundantly download gigabytes of environments into every single sample folder!
* **`default-resources`**: Uncomment these lines and define your SLURM partition/queue. For example, change `- partition=""` to `- partition="compute_node"`.
* **`singularity-args`**: By default, `-B /:/` binds your entire filesystem to the container. If your cluster restricts this, change it to bind only your specific data drives (e.g., `"-B /scratch:/scratch,-B /projects:/projects"`).

**Example of an updated `profile/config.yaml`:**
```yaml
jobs: 10
latency-wait: 60
reason: True
keep-going: True
printshellcmds: True
rerun-incomplete: True
restart-times: 2

# Slurm default resources (Uncomment and modify for your cluster)
default-resources:
  - partition="your_cluster_partition"
# - slurm="exclude=compute[011-013]"

# MUST BE ABSOLUTE PATHS
conda-prefix: /absolute/path/to/somatic_variant_caller/conda_env_singularity
singularity-prefix: /absolute/path/to/somatic_variant_caller/Singularity

# Modify this for the bindings Singularity needs on your HPC
singularity-args: "-B /:/"
```
---

## Setup and Pipeline Execution

The workflow is orchestrated using the provided helper notebook (`helper.ipynb`).
* Please follow the steps within the notebook to run the pipeline.
* Ensure the notebook is launched and run directly from the `base_env` mamba environment located at the root of your cloned repository.

### Cache Environments and Containers

On many HPC clusters, compute nodes do not have internet access, or downloading massive Singularity images on the fly can cause jobs to time out. You can force Snakemake to pre-install all Conda environments and pull all Singularity containers from your head node *without* executing the actual pipeline.

Because the pipeline's rules require a defined working directory to parse successfully, you must run this command **after** you have used the Jupyter notebook to generate your sample scripts.

Pass any **one** of the generated `.yaml` configuration files to the command to satisfy the parser:

```bash
snakemake \
  --use-conda \
  --use-singularity \
  --conda-create-envs-only \
  --profile cubi-v1
  --workflow-profile profile \
  -c 1 \
  --configfile path/to/results/view_by_pid/PATIENT_ID/somatic_mutations/SAMPLE_NAME/SAMPLE_NAME.yaml
```

*Note: You only need to do this once. Once the environments and containers are cached in your `.snakemake` folder, all subsequently submitted cluster jobs will automatically detect and use them.*

---

## Metadata Preparation

The pipeline requires a metadata dataframe to map your `.cram` or `.bam` files. If you used the ChapuyLab alignment pipeline, the notebook can automatically generate this for you.

Otherwise, prepare a TSV (Tab-Separated Values) file with the following columns:

| Column Name | Description |
| :--- | :--- |
| `BAM_FILE` | **(Required)** Absolute path to the alignment file. |
| `SAMPLE_TYPE` | **(Required)** The type of sample (e.g., "tumor", "normal", "control"). **Note:** Normal/control samples *must* contain the word "normal" or "control" in this field. Do not use underscores (`_`); use hyphens (`-`) instead. |
| `PATIENT_ID` | **(Required)** A unique identifier for the patient. Do not use underscores (`_`); use hyphens (`-`) instead. |
| `SAMPLE_NAME` | **(Required)** A unique identifier combining Patient ID and Sample Type (e.g., `PATIENT-1_tumor`). |
| `INSERT_SIZE` | **(Optional)** Absolute path to the output of GATK `CollectInsertSizeMetrics` (e.g., `..._insert_size_metrics.txt`). If omitted, the pipeline will compute this automatically, but providing it will save compute time. |

> **Safety Check:** The helper notebook enforces strict uniqueness on the `SAMPLE_NAME` column to prevent job collisions and accidental data overwriting. If any underscores (`_`) are found in `SAMPLE_TYPE` or `PATIENT_ID`, the notebook will automatically convert them to hyphens (`-`) to ensure downstream compatibility.
---

## Special Notes

> ⚠️ **Important Environment Recommendation**
> It is highly recommended that the pipeline is run in a **Singularity environment**. We have observed discrepancies in results when not using Singularity in a cluster environment.
