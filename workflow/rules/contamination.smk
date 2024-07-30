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
        normal=wrkdir / str(normal + ".pileup.table") if normal else [],
    output:
        table=temp(wrkdir / "contamination.table"),
        tumour_segments=temp(wrkdir / "tumour.segments"),
    params:
        matched=(
            "-matched " + str(wrkdir / str(normal + ".pileup.table")) if normal else ""
        ),
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
        "gatk --java-options '-Xmx{resources.mem_mb}m' CalculateContamination -I {input.tumour} {params.matched} -O {output.table} --tumor-segmentation {output.tumour_segments} &> {log}"


rule CrossCheckFingerprint:
    input:
        normal=wrkdir / "alignments" / (normal + ".bam") if normal else [],
        tumour=wrkdir / "alignments" / (tumour + ".bam"),
        haplotypemap=haplotypemap,
    params:
        LOD=LOD,
        expect_all_groups_to_match=expect_all_groups_to_match,
    output:
        crosscheckmetricsfile=wrkdir / str(output_prefix + ".crosscheck_metrics"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    log:
        logdir / "gatk/crosscheckfingerprint.log",
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' CrosscheckFingerprints -I {input.tumour} -I {input.normal} -O {output.crosscheckmetricsfile} "
        "-H {input.haplotypemap} --LOD {params.LOD} --EXPECT_ALL_GROUPS_TO_MATCH {params.expect_all_groups_to_match} --CROSSCHECK_BY SAMPLE &> {log}"
