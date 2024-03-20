rule SplitIntervals:
    input:
        genome=genome,
        target_file=target_file if target_file else [],
    params:
        target_file="-L " + target_file if target_file else "",
        scatter_count=scatter_count,
    output:
        interval_files=expand(
            wrkdir / "intervals" / "{scatter}-scattered.interval_list",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
        intervals=directory(wrkdir / "intervals"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    log:
        logdir / "gatk/splitintervals.log",
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' SplitIntervals -R {input.genome} "
        "{params.target_file} --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count {params.scatter_count} -O {output.intervals} &> {log}"


rule mutect:
    input:
        normal=wrkdir / "alignments" / (normal + ".bam"),
        tumour=wrkdir / "alignments" / (tumour + ".bam"),
        genome=genome,
        germline=germline,
        interval=wrkdir / "intervals" / "{scatter}-scattered.interval_list",
    params:
        rg_tumour=config["rg_tumour"],
        rg_normal=config["rg_normal"],
    output:
        vcf=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf"),
        stats=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf.stats"),
        f1r2=temp(wrkdir / "tmp" / "f1r2_{scatter}.tar.gz"),
        idx=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf.idx"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    log:
        logdir / "gatk/mutect_{scatter}.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' Mutect2 -R {input.genome} "
        " -I {input.normal} -I {input.tumour} "
        "-tumor {params.rg_tumour} "
        "-normal {params.rg_normal} "
        "-germline-resource {input.germline} "
        "--native-pair-hmm-threads {threads} "
        "--f1r2-tar-gz {output.f1r2} "
        "-L {input.interval} "
        "-O {output.vcf} &> {log}"


rule GatherVCFFiles:
    input:
        idx=expand(
            wrkdir / "tmp" / "unfiltered_{scatter}.vcf.idx",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
        vcf=expand(
            wrkdir / "tmp" / "unfiltered_{scatter}.vcf",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
    output:
        vcf=temp(wrkdir / "unfiltered.vcf"),
        idx=temp(wrkdir / "unfiltered.vcf.idx"),
    conda:
        "../envs/gatk.yaml"
    params:
        vcf=" ".join(
            [
                "-I " + i
                for i in expand(
                    wrkdir / "tmp" / "unfiltered_{scatter}.vcf",
                    scatter=[
                        "00" + str(i) if i > 9 else "000" + str(i)
                        for i in range(scatter_count)
                    ],
                )
            ]
        ),
    threads: 1
    log:
        logdir / "gatk/gather_vcf.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' GatherVcfs {params.vcf} -O {output.vcf} &> {log}"


rule MergeMutectStats:
    input:
        idx=expand(
            wrkdir / "tmp" / "unfiltered_{scatter}.vcf.idx",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
        stats=expand(
            wrkdir / "tmp" / "unfiltered_{scatter}.vcf.stats",
            scatter=[
                "00" + str(i) if i > 9 else "000" + str(i)
                for i in range(scatter_count)
            ],
        ),
    output:
        stats=wrkdir / "unfiltered.vcf.stats",
    conda:
        "../envs/gatk.yaml"
    params:
        stats=" ".join(
            [
                "-stats " + i
                for i in expand(
                    wrkdir / "tmp" / "unfiltered_{scatter}.vcf.stats",
                    scatter=[
                        "00" + str(i) if i > 9 else "000" + str(i)
                        for i in range(scatter_count)
                    ],
                )
            ]
        ),
    threads: 1
    log:
        logdir / "gatk/merge_stats.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' MergeMutectStats {params.stats} -O {output.stats} &> {log}"
