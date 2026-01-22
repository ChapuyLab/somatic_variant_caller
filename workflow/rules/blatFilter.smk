rule getInsertSize:
    input:
        bam=(
            tumour_bam_file
            if tumour_bam_file
            else wrkdir / "alignments" / (tumour + ".bam")
        ),
    output:
        insert_size=wrkdir / "alignments" / (tumour + "_insert-size.txt"),
        insert_size_pdf=wrkdir / "alignments" / (tumour + "_insert-size.pdf"),
    conda:
        "../envs/gatk.yaml"
    threads: 1
    resources:
        mem_mb=24000,
        runtime=24 * 60,
        nodes=1,
    log:
        logdir / "picard" / "insert_size.log",
    message:
        "Collecting insert size metrics"
    shell:
        "gatk CollectInsertSizeMetrics -I {input.bam} -O {output.insert_size} -H {output.insert_size_pdf}"  #" &> {log}"


rule gen_occ:
    input:
        genome=genome,
    params:
        file_name="11.occ",
    output:
        output_dir=directory(wrkdir / "tmp"),
        occ=temp(wrkdir / "tmp" / "11.occ"),
        blat_temp_out=temp(wrkdir / "tmp" / "genome.psl"),
    threads: 1
    log:
        logdir / "blat" / "Gen_Occ.log",
    resources:
        mem_mb=4000,
        runtime=20,
        nodes=1,
    conda:
        "../envs/blatfilter.yaml"
    shell:
        "cd {output.output_dir}; "
        "blat "
        "{input.genome} "
        "{input.genome} "
        "-makeOoc={params.file_name} "
        "-t=dna -q=dna "
        "{output.blat_temp_out} "


rule gen_2bit:
    input:
        genome=genome,
    output:
        temp(wrkdir / "tmp" / "genome.2bit"),
    threads: 1
    log:
        logdir / "tmp" / "Gen_2Bit_genome.log",
    resources:
        mem_mb=4000,
        runtime=20,
        nodes=1,
    conda:
        "../envs/blatfilter.yaml"
    shell:
        "faToTwoBit "
        "{input.genome} "
        "{output} "


checkpoint scatter_maf:
    input:
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
    params:
        max_mut=10000,
    output:
        maf_tmp=directory(wrkdir / "scatter_tmp_{sample}"),
    threads: 1
    log:
        logdir / "blat" / "scatter_{sample}.log",
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
    run:
        import numpy as np
        import pandas as pd
        from pathlib import Path
        import os

        set_input = set()
        line_count = 0
        skipline = 0
        otherinfoIndex = ""
        with open(input.maf) as handle:
            for line in handle:
                if "#" in line or "Hugo_Symbol" in line:
                    if "Hugo_Symbol" in line:
                        otherinfoIndex = line.strip().split("\t").index("Otherinfo2")
                    skipline += 1
                    continue
                else:
                    line = line.strip().split("\t")
                    set_input.add(line[4] + ":" + line[otherinfoIndex])
                    line_count += 1
                    # Calculate # of Chunks
        if line_count > params.max_mut:
            scatter_count = line_count // params.max_mut
        else:
            scatter_count = 1
        chunksize = np.array_split(np.array(range(0, line_count)), scatter_count)
        # print(chunksize)
        # print(line_count)
        ## Writeout maf in smaller chunks

        scatter_index = 0
        start = 0
        header = None
        set_output = set()
        outdir = Path(output.maf_tmp)
        os.makedirs(outdir, exist_ok=True)
        # for df_maf in pd.read_csv(input.maf, sep="\t", low_memory=False, chunksize=len(chunksize[scatter_index]), comment='#'):

        for size in chunksize:
            df_maf = pd.read_csv(
                input.maf,
                skiprows=start,
                nrows=len(size),
                comment="#",
                header=0 if start == 0 else None,
                names=header if start != 0 else None,
                sep="\t",
            )
            if start == 0:
                header = df_maf.columns
                start += len(size) + 1
            else:
                start += len(size)
            set_output = set_output.union(
                set(
                    df_maf["Chromosome"].astype(str)
                    + ":"
                    + df_maf["Otherinfo2"].astype(str)
                )
            )

            df_maf.Start_Position = df_maf.Start_Position.astype(float).astype(int)
            df_maf.End_Position = df_maf.End_Position.astype(float).astype(int)

            df_out = df_maf.to_csv(
                outdir / (wildcards.sample + "_" + str(scatter_index) + ".maf"),
                index=False,
                sep="\t",
            )
            scatter_index += 1
        if len(set_input - set_output) > 0:
            print(line_count)
            print(params.max_mut)
            print(scatter_count)
            print(set_input - set_output)
            raise ValueError("Error")


rule BlatFilter:
    input:
        maf=lambda wildcards: wrkdir
        / ("scatter_tmp_" + wildcards.sample)
        / (wildcards.sample + "_" + wildcards.scatter + ".maf"),
        # f"scatter_tmp_{wildcards.sample}" / f"{wildcards.sample}_{wildcards.scatter}.maf",
        bam=(
            tumour_bam_file
            if tumour_bam_file
            else wrkdir / "alignments" / (tumour + ".bam")
        ),
        bai=(
            tumour_bam_file
            if tumour_bam_file
            else wrkdir / "alignments" / (tumour + ".bam.bai")
        ),
        insert_size=(
            wrkdir / "alignments" / (tumour + "_insert-size.txt")
            if insertSize is None
            else insertSize
        ),
        database=wrkdir / "tmp" / "genome.2bit",
        occ=wrkdir / "tmp" / "11.occ",
        genome=genome,
    params:
        blat_binary="blat",
        output_prefix=lambda wildcards: str(
            wrkdir / "tmp_filter" / (wildcards.sample + "_" + wildcards.scatter)
        ),
        output_dir=str(wrkdir / "tmp_filter"),
        stepper=stepper,
    output:
        mate_query=temp(wrkdir / "tmp_filter" / "{sample}_{scatter}_mate_query.fa"),
        psl=temp(wrkdir / "tmp_filter" / "{sample}_{scatter}_output.psl"),
        query=temp(wrkdir / "tmp_filter" / "{sample}_{scatter}_query.fa"),
        maf_all=temp(wrkdir / "tmp_filter" / "{sample}_{scatter}.blat.all.maf"),
        maf_rejected=temp(
            wrkdir / "tmp_filter" / "{sample}_{scatter}.blat.rejected.maf"
        ),
        maf=temp(wrkdir / "tmp_filter" / "{sample}_{scatter}.blat.passed.maf"),
    threads: 1
    conda:
        "../envs/blatfilter.yaml"
    log:
        logdir / "blat" / "{sample}_{scatter}_blat_filter.log",
    resources:
        mem_mb=4000,
        runtime=60 * 24 * 1,
        nodes=1,
    script:
        "../scripts/blatfilter/blat_filter.py"


def getMafs(wildcards):
    scatter_dir = Path(checkpoints.scatter_maf.get(**wildcards).output[0])
    scatter = glob_wildcards(
        scatter_dir / (wildcards.sample + "_{scatter}.maf")
    ).scatter
    files = list()
    filter_dir = Path(rules.BlatFilter.params.output_dir)
    for i in scatter:
        files.append(
            filter_dir
            / str(wildcards.sample + "_" + str(i) + ".blat." + wildcards.group + ".maf")
        )
    return files


rule gather_maf:
    input:
        mafs=getMafs,
    params:
        chunksize=10000,
    output:
        maf=wrkdir / "blat" / "{sample}_{group}.maf",
    threads: 1
    log:
        logdir / "blat" / "{sample}_{group}_gather_maf.log",
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
    run:
        import pandas as pd

        header = False
        with open(output.maf, "w") as handle:
            for i in input.mafs:
                line_count = 0
                with open(i) as handle_in:
                    for line in handle_in:
                        line_count += 1
                        if line_count > 1:
                            break
                if line_count == 1:
                    continue
                if not header:
                    with open(i) as handle_in:
                        for line in handle_in:
                            handle.write(line)
                            header = True
                            break
                for df_maf in pd.read_csv(
                    i,
                    sep="\t",
                    skiprows=1,
                    header=None,
                    low_memory=False,
                    chunksize=params.chunksize,
                ):
                    df_out = df_maf.to_csv(index=False, header=False, sep="\t")
                    handle.write(df_out)


rule FilterVaf:
    input:
        maf=wrkdir / "blat" / str(output_prefix + "_rejected.maf"),
        vcf=wrkdir / str(output_prefix + "{ifrescue}.vcf"),
    params:
        chunksize=10000,
    output:
        vcf=wrkdir / str(output_prefix + "{ifrescue}.blat.vcf"),
    threads: 1
    log:
        logdir / "blat" / "{ifrescue}_gather_maf.log",
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
    run:
        set_reject = set()
        with open(input.maf) as handle:
            for line in handle:
                if "#" in line or "Hugo_Symbol" in line:
                    otherinfoIndex = line.strip().split("\t").index("Otherinfo1")
                    otherinfo2Index = line.strip().split("\t").index("Otherinfo2")
                    continue
                else:
                    line = line.strip().split("\t")
                    if line[-1] == "REJECT":
                        set_reject.add(
                            line[otherinfoIndex] + ":" + line[otherinfo2Index]
                        )
        print("Number of mutations to be filtered out ", len(set_reject))
        with open(input.vcf) as handle, open(output.vcf, "w") as handle_out:
            flag = True
            for line in handle:
                line = line.strip().split("\t")
                if "#" == line[0][0]:
                    if "##FILTER" in line and flag:
                        handle_out.write(
                            '##FILTER=<ID=blatReject,Description="Mutation does not meet criteria of a blat Filter">'
                        )
                        flag = False
                    if "##filtering_status" in line:
                        handle_out.write(
                            "##filtering_status=These calls have been further filtered using a blat filter"
                        )
                elif line[0] + ":" + line[1] in set_reject:
                    line[6] = "blatReject"
                handle_out.write("\t".join(line) + "\n")
