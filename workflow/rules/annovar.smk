if run_annovar:

    rule convert2annovar:
        input:
            vcf=wrkdir / str(output_prefix + "_filtered.vcf"),
        output:
            anno_norm=(
                wrkdir / "AV_input" / ("AV." + rg_normal + ".avinput")
                if normal
                else {}
            ),
            anno_tum=wrkdir / "AV_input" / ("AV." + rg_tumour + ".avinput"),
            annovar=directory(wrkdir / "AV_input"),
        params:
            annovar=wrkdir / "AV_input" / "AV",
        conda:
            "../envs/annovar.yaml"
        threads: 1
        log:
            logdir / "annovar/convert2annovar.log",
        resources:
            mem_mb=1000,
            runtime=72 * 60,
            nodes=1,
        shell:
            (
                "perl "
                + str(path_to_annovar / "convert2annovar.pl")
                + " -filter pass -format vcf4 -allsample {input.vcf} -outfile {params.annovar} &> {log}"
            )

    rule annotate:
        input:
            annovar=wrkdir / "AV_input" / ("AV." + rg_tumour + ".avinput"),
        output:
            annotated=wrkdir
            / str(output_prefix + "_filtered.avinput.hg19_multianno.csv"),
        params:
            annotated=wrkdir / str(output_prefix + "_filtered.avinput"),  # this is need as annovar will add the suffix to the output file
            genome_ver="hg19",
            protocol="refGene,cytoBand,avsnp150,dbnsfp42c",
            operation="g,r,f,f",
        conda:
            "../envs/annovar.yaml"
        threads: 1
        log:
            logdir / "annovar/annotate.log",
        resources:
            mem_mb=1000,
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
                + " -operation {params.operation} -nastring . -csvout -polish &> {log}"
            )
