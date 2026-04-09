
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
        table=wrkdir / "read-orientation-model.tar.gz",
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
        mem_mb=20000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' LearnReadOrientationModel {params.f1r2} -O {output.table} &> {log}"


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
        vcf=wrkdir / str(output_prefix + "_filtered.vcf"),
        idx=wrkdir / str(output_prefix + "_filtered.vcf.idx"),
        stats=wrkdir / str(output_prefix + "_filtered.vcf.filteringStats.tsv"),
    params:
        filter_mutect_extra=(
            config["filter_mutect_extra"] if "filter_mutect_extra" in config else ""
        ),
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
        "--ob-priors {input.read_orientation_model} {params.filter_mutect_extra} &> {log}"


# rule vcf2maf:
#     input:
#         vcf=wrkdir / str(output_prefix + "_filtered.vcf")


rule convert2callstats:
    input:
        vcf=wrkdir / str(output_prefix + "_filtered.vcf"),
    output:
        callstats=temp(wrkdir / str(output_prefix + "_callstats.txt")),
    conda:
        "../envs/detin.yaml"
    threads: 1
    log:
        logdir / "detin/callstats_convertor.log",
    resources:
        mem_mb=50000,
        runtime=72 * 60,
        nodes=1,
    script:
        "../scripts/vcf2callstats.py"


rule detin:
    input:
        mutation_data_path=wrkdir / str(output_prefix + "_callstats.txt"),
        indel_data_path=wrkdir / str(output_prefix + "_filtered.vcf"),
        exac_data_path=pickle_high_af,
    params:
        output_name=output_prefix,
        indel_data_type="mutect2",
    output:
        snvs=wrkdir / "detin" / str(output_prefix + ".deTiN_SSNVs.txt"),
        indels=wrkdir / "detin" / str(output_prefix + ".deTiN_indels.txt"),
        output_dir=directory(wrkdir / "detin"),
    conda:
        "../envs/detin.yaml"
    threads: 1
    log:
        logdir / "detin/detin.log",
    resources:
        mem_mb=50000,
        runtime=72 * 60,
        nodes=1,
    script:
        "../scripts/detin_scripts/deTiN.py"
        # "--mutation_data_path {input.callstats} "
        # "--exac_data_path {input.pickle_high_af} "
        # "--indel_data_path {input.vcf} "
        # "--output_name {params.output_name} "
        # "--output_dir {params.output_path} "
        # "--indel_data_type mutect2"


rule rescueMutations:
    input:
        snvs=wrkdir / "detin" / str(output_prefix + ".deTiN_SSNVs.txt"),
        indels=wrkdir / "detin" / str(output_prefix + ".deTiN_indels.txt"),
        vcf=wrkdir / str(output_prefix + "_filtered.vcf"),
    output:
        vcf=wrkdir / str(output_prefix + "_filtered.rescue.vcf"),
    threads: 1
    conda:
        "../envs/detin.yaml"
    log:
        logdir / "detin/rescue.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    script:
        "../scripts/detin_scripts/flip2pass.py"
