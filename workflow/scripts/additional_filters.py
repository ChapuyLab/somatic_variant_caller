import os
import pandas as pd
import io

def read_vcf(path):
    lines=list()
    header=list()
    with open(path, 'r') as f:
        for l in f:
            if l.startswith('##'):
                header.append(l.strip())
            else:
                lines.append(l)
    return header, pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t'
    )

def write_vcf(header, dataframe, path):
    # print(header)
    with open(path, 'w') as handle:
        handle.write('\n'.join(header))
        handle.write('\n')
    with open(path, 'a') as handle:
        dataframe.to_csv(handle, index=False, sep='\t')


def is_empty(info, key):
    if key in info:
        if info[key] == "." or info[key] == "0":
            return True
        else:
            return False
    return True

def annotate_confidence(INFO):
    confidence=10
    in1KG = False
    indbSNP = False
    precious = False
    is_repeat = False
    is_STR = False
    is_weird = False
    # print(INFO)
    INFO=dict(item.split("=") for item in INFO.split(";") if '=' in item)
    # print(INFO)
    if 'REPET' in INFO:
        # print(INFO)
        if ('Simple_repeat' in INFO['REPET'] or 'Low_' in INFO['REPET'] or 'Satellite' in INFO['REPET']):
            is_repeat=True
            confidence-=2
        elif not is_empty(INFO, 'REPET'):
            confidence-=1
    if 'STREP' in INFO:
        if not is_empty(INFO, 'STREP'):
            is_STR=True
            if not is_repeat:
                confidence-=2
    if (not is_empty(INFO, 'EXCLU')) or (not is_empty(INFO, 'BLACK')):
        confidence=-3
        is_weird=1
    if not is_empty(INFO, 'HSDEP'):
        confidence-=3
        is_weird=1
    if 'MAPAB' in INFO:
        if is_empty(INFO, 'MAPAB'):
            confidence -=5
        else:
            mapab=float(INFO['MAPAB'])
            if  mapab < 0.5:
                confidence -=1
                is_weird=True
                if mapab < 0.4:
                    confidence-=1
                    if mapab <0.25:
                        confidence=-1
                        if mapab<0.1:
                            confidence-=2
                            if mapab<0.005:
                                confidence-=3
    else:
        confidence=-5
	# if others have found the SNP already, it may be interesting despite low score - but only if it's not a weird region.
	# if a position gets up from score 4 to 6 by giving +2 for presence in dbSNP, it's very likely an artefact reported to dbSNP
	# if (($in1KG || $indbSNP) && ! $is_weird)
	# {
	# 	$confidence++;
	# }
    return_tuple=['', '']
    if confidence<1:
        INFO['QC_CONF']="1"
    elif confidence>10:
        INFO['QC_CONF']="10"
    else:
        INFO['QC_CONF']=str(confidence)

    return_tuple[1]= ";".join(f"{key}={value}" for key, value in INFO.items())

    if confidence >= 8:
        return_tuple[0]='PASS'
    else:
        return_tuple[0]='QC_ISSUE'
    return return_tuple

import argparse

parser = argparse.ArgumentParser(
                    prog='additonal_filters.py',
                    description='Assigns a confidence measure to each variant detected by mutect2',
                    )


parser.add_argument('filename')

args = parser.parse_args()
header, vcf_df = read_vcf(args.filename)
vcf_df[['FILTER', 'INFO']]=vcf_df.INFO.apply(annotate_confidence).apply(pd.Series)

write_vcf(header, vcf_df, args.filename)

print('Total Variants before:' ,vcf_df.shape[0])
print('Total Variants after:' ,vcf_df[vcf_df.FILTER!='QC_ISSUE'].shape[0])
