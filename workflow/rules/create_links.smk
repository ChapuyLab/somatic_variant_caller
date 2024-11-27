rule create_links:
    input:
        metadata=metadata,
    params:
        alignment_dir=wrkdir / "alignments",
    output:
        expand(wrkdir / "alignments" / "{sample}.bam", sample=samples),  # these need to be temp but making them temp creates a problem at the checkpoint hence letting them be as they are just links NEED TO FIX
        expand(wrkdir / "alignments" / "{sample}.bam.bai", sample=samples),
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
    threads: 1
    run:
        import pandas as pd

        metadata = pd.read_csv(input.metadata)
        metadata = metadata[
            (metadata["SAMPLE_TYPE"].isin(samples))
            & (metadata["PATIENT_ID"].astype(str) == str(config["pid"]))
        ]
        for index, row in metadata.iterrows():
            fastq_file = Path(row["BAM_FILE"])
            output_file = params.alignment_dir / (row["SAMPLE_TYPE"] + ".bam")
            # os.makedirs(params.fastq_dir, exist_ok=True)
            if os.path.exists(output_file):
                os.remove(output_file)
            os.symlink(fastq_file, output_file)
            if os.path.exists(str(output_file) + ".bai"):
                os.remove(str(output_file) + ".bai")
            os.symlink(str(fastq_file) + ".bai", str(output_file) + ".bai")
