rule pileupSummaries:
    input:
        bam=wrkdir / "alignments" / "{sample}.bam",
        variant=config["common_biallelic"],
    output:
        table=temp(wrkdir / "{sample}.pileup.table"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    log:
        logdir / "gatk/{sample}_pileup.log",
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' GetPileupSummaries -I {input.bam} -V {input.variant} -L {input.variant} -O {output.table} &> {log}"


rule CalculateContamination:
    input:
        tumour=wrkdir / str(tumour + ".pileup.table"),
        normal=wrkdir / str(normal + ".pileup.table"),
    output:
        table=temp(wrkdir / "contamination.table"),
        tumour_segments=temp(wrkdir / "tumour.segments"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    log:
        logdir / "gatk/calculatecontamination.log",
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' CalculateContamination -I {input.tumour} -matched {input.normal} -O {output.table} --tumor-segmentation {output.tumour_segments} &> {log}"
