rule annovar2maf:
    input:
        annotated=wrkdir / str(output_prefix + "{ifrescue}.avinput.hg19_multianno.txt"),
    output:
        maf=wrkdir / str(output_prefix + "{ifrescue}.avinput.hg19_multianno.maf"),
    params:
        tsb=output_prefix,
        build=genome_build,
        protocol=geneDatabase,  # replace this later
    conda:
        "../envs/detin.yaml"
    threads: 1
    log:
        logdir / "createMaf/annovar2maf{ifrescue}.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    script:
        "../scripts/annovar2maf.py"


rule vcf2maf:
    input:
        vcf=wrkdir / str(output_prefix + "{ifrescue}.vcf"),
    output:
        maf=wrkdir / str(output_prefix + "{ifrescue}.maf"),
    params:
        inhibit_vep="--inhibit-vep ",
        genome=genome,
        normal_id=" --normal-id " + rg_normal if normal else "",
        tumor_id=" --tumor-id " + rg_tumour,
        protocol="refGene",  # replace this later
    conda:
        "../envs/detin.yaml"
    threads: 1
    log:
        logdir / "createMaf/vcf2maf{ifrescue}.log",
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    shell:
        "vcf2maf.pl "
        "{params.inhibit_vep}"
        "--input-vcf {input.vcf} "
        "--output-maf {output.maf} "
        "--ref-fasta {params.genome} "
        "{params.tumour_id} {params.normal_id} &> {log}"


rule addExtraFields:
    input:
        maf=(
            wrkdir / str(output_prefix + "_filtered.rescue.avinput.hg19_multianno.maf")
            if normal
            else wrkdir / str(output_prefix + "_filtered.avinput.hg19_multianno.maf")
        ),
    output:
        maf=(
            wrkdir
            / str(
                output_prefix
                + "_filtered.rescue.avinput.hg19_multianno.extended_info.maf"
            )
            if normal
            else wrkdir
            / str(output_prefix + "_filtered.avinput.hg19_multianno.extended_info.maf")
        ),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=24 * 60,
        nodes=1,
    log:
        logdir / "annovar" / "extended_annotations.log",
    message:
        "Adding VAF, ref count information to maf"
    run:
        with open(input.maf) as handle, open(output.maf, "w") as handle_out:
            for line in handle:
                if "#" in line:
                    line = line.strip()
                elif "Hugo_Symbol" in line:
                    line = line.strip().split("\t")
                    line += ["t_ref_count", "t_alt_count", "AF"]
                else:
                    # print(line)
                    line = line.strip().split("\t")
                    # print(len(line))
                    info_col = line[129].split(":")
                    line += [
                        info_col[1].split(",")[0],
                        info_col[1].split(",")[1],
                        info_col[2],
                    ]
                handle_out.write("\t".join(line) + "\n")
