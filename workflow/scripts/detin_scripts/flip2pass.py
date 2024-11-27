flip_to_pass=set()

with open(snakemake.input['snvs']) as handle:
    for line in handle:
        line=line.strip().split()
        if 'contig' == line[0]:
            continue
        else:
            if line[7] == 'KEEP':
                flip_to_pass.add(str(line[0])+':'+str(line[1]))

with open(snakemake.input['indels']) as handle:
    for line in handle:
        line=line.strip().split()
        if 'contig' == line[0]:
            continue
        else:
            if line[6] == 'PASS' and line[6] != line[11]:
                flip_to_pass.add(str(line[0])+':'+str(line[1]))

with open(snakemake.input['vcf']) as handle, open(snakemake.output['vcf'], 'w') as handle_out:
    for line in handle:
        if '##' in line:
            pass
        elif '#CHROM' in line:
            handle_out.write('##detinRescueStatus=These calls have been rescued using detin SSNV only mode\n')
        else:
            line=line.strip().split('\t')
            if line[0]+':'+line[1] in flip_to_pass:
                line[6] = 'PASS'
            line='\t'.join(line)+'\n'
        handle_out.write(line)
