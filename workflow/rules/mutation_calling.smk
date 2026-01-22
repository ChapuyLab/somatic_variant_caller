rule SplitIntervals:
    input:
        genome=genome,
        target_file=target_file if target_file else [],
    params:
        target_file=(
            "-L " + target_file
            if target_file
            else "-L /dh-projects/ag-ishaque/analysis/sahays/pipelines/mutect2_pipeline/reference/wgs_interval_file.bed"
        ),
        scatter_count=scatter_count,
        subdivision_mode=(
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
            if target_file
            else "INTERVAL_COUNT"
        ),
    output:
        interval_files=temp(
            expand(
                wrkdir / "intervals" / "{scatter}-scattered.interval_list",
                scatter=[
                    "00" + str(i) if i > 9 else "000" + str(i)
                    for i in range(scatter_count)
                ],
            )
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
        "{params.target_file} --subdivision-mode {params.subdivision_mode} "
        "--scatter-count {params.scatter_count} -O {output.intervals}"
        #" &> {log}"


def getInterval(scatter):
    if len(str(int(scatter) + 1)) == 1:
        return "0" + str(int(scatter) + 1)
    else:
        return str(int(scatter) + 1)


def getIntervalFile(wildcards):
    intervalFile = (
        interval_folder
        / (getInterval(wildcards.scatter) + "_of_80")
        / "scattered.interval_list"
    )
    return intervalFile


rule mutect:
    input:
        normal=(
            normal_bam_file
            if normal_bam_file
            else wrkdir / "alignments" / (normal + ".bam") if normal else []
        ),
        tumour=(
            tumour_bam_file
            if tumour_bam_file
            else wrkdir / "alignments" / (tumour + ".bam")
        ),
        genome=genome,
        germline=germline,
        interval=(
            getIntervalFile
            if interval_folder
            else wrkdir / "intervals" / "{scatter}-scattered.interval_list"
        ),
        pon=pon if pon else [],
    params:
        normal=(
            "-I "
            + (
                normal_bam_file
                if normal_bam_file
                else wrkdir / "alignments" / (normal + ".bam")
            )
            if normal
            else []
        ),
        rg_tumour=config["rg_tumour"],
        rg_normal="-normal " + config["rg_normal"] if normal else "",
        pon="-pon " + pon + " --genotype-pon-sites true " if pon else "",
    output:
        vcf=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf"),
        stats=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf.stats"),
        f1r2=temp(wrkdir / "tmp" / "f1r2_{scatter}.tar.gz"),
        idx=temp(wrkdir / "tmp" / "unfiltered_{scatter}.vcf.idx"),
        # tmpdir=temp(directory(wrkdir / "tmp")),
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
        "-I {input.tumour} "
        "{params.normal} "
        "{params.rg_normal} "
        "--germline-resource {input.germline} "
        "--native-pair-hmm-threads {threads} "
        "--genotype-germline-sites true "
        "--f1r2-tar-gz {output.f1r2} "
        "{params.pon} "
        "-L {input.interval} "
        "-O {output.vcf}"
        #" &> {log}"


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
        vcf=wrkdir / "unfiltered.vcf",
        idx=wrkdir / "unfiltered.vcf.idx",
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
        "gatk --java-options '-Xmx{resources.mem_mb}m' GatherVcfs {params.vcf} -O {output.vcf}"
        #" &> {log}"


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
        "gatk --java-options '-Xmx{resources.mem_mb}m' MergeMutectStats {params.stats} -O {output.stats}"
        #" &> {log}"
