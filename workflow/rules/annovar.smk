rule convert2annovar:
    input:
        vcf=wrkdir / str(output_prefix + "{ifrescue}.vcf"),
    output:
        anno_tum=temp(
            wrkdir / "AV_input{ifrescue}" / ("AV{ifrescue}." + rg_tumour + ".avinput")
        ),
        annovar=directory(wrkdir / "AV_input{ifrescue}"),
    params:
        annovar=lambda wildcards: wrkdir
        / str("AV_input" + wildcards.ifrescue)
        / str("AV" + wildcards.ifrescue),
    conda:
        "../envs/annovar.yaml"
    threads: 1
    log:
        logdir / "annovar/convert2annovar{ifrescue}.log",
    resources:
        mem_mb=1000,
        runtime=72 * 60,
        nodes=1,
    shell:
        (
            "perl "
            + str(path_to_annovar / "convert2annovar.pl")
            + " -filter pass -includeinfo -format vcf4 -allsample {input.vcf} -outfile {params.annovar} &> {log}"
        )


rule annotate:
    input:
        annovar=wrkdir
        / "AV_input{ifrescue}"
        / ("AV{ifrescue}." + rg_tumour + ".avinput"),
    output:
        annotated=temp(
            wrkdir / str(output_prefix + "{ifrescue}.avinput.hg19_multianno.txt")
        ),
    params:
        annotated=lambda wildcards: wrkdir
        / str(output_prefix + wildcards.ifrescue + ".avinput"),
        # this is need as annovar will add the suffix to the output file
        genome_ver="hg19",
        protocol="refGene,cytoBand,avsnp150,dbnsfp42c",
        operation="g,r,f,f",
    conda:
        "../envs/annovar.yaml"
    threads: 1
    log:
        logdir / "annovar/annotate{ifrescue}.log",
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
    shell:
        (
            "perl "
            + str(path_to_annovar / "table_annovar.pl")
            + " {input.annovar} "
            + str(annovar_db)
            + " "
            + "-buildver {params.genome_ver} -outfile {params.annotated} -remove -protocol {params.protocol}"
            + " -operation {params.operation} --otherinfo -nastring . -polish &> {log}"
        )
