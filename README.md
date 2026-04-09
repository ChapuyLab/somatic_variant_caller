# ChapuyLab Somatic Variant Caller


This repository contains a Snakemake-based pipeline for running somatic variant calling using GATK Mutect2. It is designed to process individual samples (both tumor-normal paired and tumor-only) efficiently on HPC clusters.

The setup and job orchestration are managed via a helper Jupyter Notebook, which generates sample-specific configurations and execution scripts.

## Features

Sample Configurations: Supports matched Tumor-Normal pairs, Tumor-Only runs, and "Fake Normal" fallbacks.

Cluster Ready: Pre-configured for SLURM execution, but easily adaptable to other scheduling engines via Snakemake profiles.

Annotation & Filtering: Integrates ANNOVAR for variant annotation and deTiN for estimating tumor-in-normal contamination.

Containerized: Leverages Conda/Mamba and Singularity for reproducible environments.

Tested Build: The pipeline has been extensively tested on the hg19 genome build.


## Prerequisites and Required Data

Before running the setup notebook, you must download and configure several reference databases and tools.

1. Reference Files (GATK Best Practices)
Download the following from the Google Cloud GATK Best Practices bucket:

Germline resource (af-only-gnomad.raw.sites.vcf)

Common biallelic variants (small_exac_common_3.vcf)

Panel of Normals (pon.vcf)

Reference Genome (hs37d5_PhiX.fa or the standard Broad hg19 build)

Haplotype map (hg19_nochr.map) from https://github.com/naumanjaved/fingerprint_maps/tree/master/map_files


2. Annotation & Filtering Resources

ANNOVAR: You must download the ANNOVAR software and its human database (humandb).

deTiN Data: Requires high AF ExAC pickle file (exac.pickle_high_af).

Target Baits: An Agilent SureSelectv2 bed file is used by default (Agilent_SureSelectv2_baits.bed).

## Setup and Pipeline Execution
The workflow is orchestrated using the provided helper notebook (helper.ipynb). Please follow the steps in there to run the pipeline Which can be run in the base_env mamba environment

## Metadata Preparation
The pipeline requires a metadata dataframe to map your .cram or .bam files. If you used the ChapuyLab alignment pipeline, the notebook can automatically generate this.

Otherwise, prepare a TSV file with the following exact columns:

BAM_FILE: Absolute path to the alignment file.

SAMPLE_TYPE: The type of sample (e.g., "tumor", "normal", "control"). Note: Control/normal samples must contain the word "normal" or "control" in this field. Do not use underscores (_).

PATIENT_ID: A unique identifier for the patient. Do not use underscores (_).

SAMPLE_NAME: A unique identifier combining Patient ID and Sample Type (e.g., PATIENT-1_tumor).

Safety Check: The notebook enforces strict uniqueness on the SAMPLE_NAME column to prevent job collisions and accidental data overwriting.

## Special Notes

It is highly recommned that the pipeline is run in singularity environment, we have observed discrepancies when not using singularity in cluster environment
