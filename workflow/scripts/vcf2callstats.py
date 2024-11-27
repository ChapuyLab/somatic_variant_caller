def is_member(a, b):
    # based on the matlab is_member function
    # code from stackoverflow user Joe Kington
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, np.nan) for itm in a]

def formatMutect2Output(vcf):
    indel_type='mutect2'
    content = []
    if vcf[-2:] == 'gz':
        with gzip.open(vcf, 'r') as f:
            content = f.readlines()
    else:
        with open(vcf) as f:
            content = f.readlines()

    cols_type = {0: str}
    for line in content:
        if line[0] == '#' and line[1] != '#':
            headerline = line.split('\t')
            break
    indel_table = pd.read_csv(vcf, sep='\t', comment='#', header=None, low_memory=False, dtype=cols_type)
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
    normal_sample = 'normal'
    tumor_sample = 'tumor'
    for line in content:
        if line[0:15] == '##normal_sample':
            normal_sample = line.split('=')[1][0:-1]
        if line[0:14] == '##tumor_sample':
            tumor_sample = line.split('=')[1][0:-1]
    if tumor_sample == 'tumor' and normal_sample == 'normal':
        indel_table.rename(
            columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'ref_allele', 4: 'alt_allele', 5: 'QUAL', 7: 'INFO', 8: 'format',
                     6: 'filter', 9: 'tumor', 10: 'normal'},
            inplace=True)
    else:
        if tumor_sample == headerline[9]:
            indel_table.rename(
                    columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'ref_allele', 4: 'alt_allele', 5: 'QUAL', 7: 'INFO', 8: 'format',
                     6: 'filter', 9: 'tumor', 10: 'normal'},
                    inplace=True)
        elif tumor_sample == headerline[10][0:-1]:
            indel_table.rename(
                columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'ref_allele', 4: 'alt_allele', 5: 'QUAL', 7: 'INFO', 8: 'format',
                         6: 'filter', 9: 'normal', 10: 'tumor'},
                inplace=True)
        else:
            raise Exception('failed to read MuTect 2 indels VCF')
    counts_format = indel_table['format'][0].split(':')
    depth_ix = counts_format.index('AD')
    indel_table = indel_table[np.isfinite(is_member(indel_table['filter'], ['PASS', 'normal_artifact']))]
    indel_table.reset_index(inplace=True, drop=True)


    # parsing format line and file to determine required alt and ref columns
    # we use "tier 1" read counts for varaints
    n_depth = np.zeros([len(indel_table), 1])
    n_alt_count = np.zeros([len(indel_table), 1])
    n_ref_count = np.zeros([len(indel_table), 1])

    t_depth = np.zeros([len(indel_table), 1])
    t_alt_count = np.zeros([len(indel_table), 1])
    t_ref_count = np.zeros([len(indel_table), 1])
    t_lod=np.zeros([len(indel_table), 1])

    for index, row in indel_table.iterrows():
        spl_n = row['normal'].split(':')
        spl_t = row['tumor'].split(':')
        n_alt_count[index] = int(spl_n[depth_ix].split(',')[1])
        n_ref_count[index] = int(spl_n[depth_ix].split(',')[0])
        n_depth[index] = n_alt_count[index]+n_ref_count[index]
        t_alt_count[index] = int(spl_t[depth_ix].split(',')[1])
        t_ref_count[index] = int(spl_t[depth_ix].split(',')[0])
        t_depth[index] = t_alt_count[index] + t_ref_count[index]
        info_dict=dict()
        for info in row.INFO.split(';'):
            if '=' in info:
                info_dict[info.split('=')[0]]=info.split('=')[1]

        t_lod[index]=info_dict['TLOD']

    # print(normal_sample)


    indel_table['t_depth'] = t_alt_count + t_ref_count
    indel_table['t_alt_count'] = t_alt_count
    indel_table['t_ref_count'] = t_ref_count

    indel_table['n_depth'] = n_alt_count + n_ref_count
    indel_table['n_alt_count'] = n_alt_count
    indel_table['n_ref_count'] = n_ref_count
    indel_table.reset_index(inplace=True, drop=True)
    indel_table['failure_reasons']=indel_table['filter'].apply(lambda x: x if x!='PASS' else '')
    indel_table['judgement']=indel_table['filter'].apply(lambda x: 'KEEP' if x=='PASS' else 'REJECT')
    indel_table['tumor_name']=tumor_sample
    indel_table['normal_name']=normal_sample
    indel_table['t_lod']=t_lod


    # annotate with acs data
    return indel_table

import argparse
import pandas as pd
import numpy as np
import gzip


# parser = argparse.ArgumentParser(
#     prog='vcf2callstats',
#     description='',
#     epilog=''
# )

# parser.add_argument('input_vcf')           # positional argument
# parser.add_argument('output_vcf')      # option that takes a value

# args=parser.parse_args()

# print(args.input_vcf)
# print(args.output_vcf)


formatMutect2Output(snakemake.input['vcf'])[
    [
        'contig','position',
        'n_alt_count',
        'n_ref_count',
        't_alt_count',
        't_ref_count',
        'failure_reasons','judgement', 't_lod', 'tumor_name', 'normal_name', 'ref_allele', 'alt_allele'
    ]
].to_csv(snakemake.output['callstats'], sep='\t' , index=False)
