
rule learnReadOrientationModel:
    input:
        f1r2=expand(
            wrkdir / "tmp" / "f1r2_{scatter}.tar.gz",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
    output:
        table=temp(wrkdir / "read-orientation-model.tar.gz"),
    conda:
        "../envs/gatk.yaml"
    params:
        f1r2=" ".join(
            [
                "-I " + i
                for i in expand(
                    wrkdir / "tmp" / "f1r2_{scatter}.tar.gz",
                    scatter=[
                        "00" + str(i) if i > 9 else "000" + str(i)
                        for i in range(scatter_count)
                    ],
                )
            ]
        ),
    threads: 4
    log:
        logdir / "gatk/readorientation.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        """
        gatk --java-options '-Xmx{resources.mem_mb}m' LearnReadOrientationModel {params.f1r2} -O {output.table} &> {log}
        """


rule filter_vcf:
    input:
        vcf=wrkdir / "unfiltered.vcf",
        stats=wrkdir / "unfiltered.vcf.stats",
        idx=wrkdir / "unfiltered.vcf.idx",
        contamination_table=wrkdir / "contamination.table",
        tumour_segments=wrkdir / "tumour.segments",
        read_orientation_model=wrkdir / "read-orientation-model.tar.gz",
        genome=genome,
    output:
        vcf=wrkdir / "filtered.vcf",
    conda:
        "../envs/gatk.yaml"
    threads: 4
    log:
        logdir / "gatk/filter_vcf.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' FilterMutectCalls -V {input.vcf} -O {output.vcf} "
        "-R {input.genome} "
        "--contamination-table {input.contamination_table} "
        "--tumor-segmentation {input.tumour_segments} "
        "--ob-priors {input.read_orientation_model} &> {log}"
