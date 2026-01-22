import argparse
from annovar2maf import read_annovar_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert annovar and bcftools-csq annotations to MAF", prog="annovar2maf")
    # parser.add_argument("-t", "--tsb", help="Sample name. Default parses from the file name")
    parser.add_argument("-i", "--input", help="Path to input multianno file", required =True)
    parser.add_argument("-o", "--output", help="Path output maf file", required =True)
    parser.add_argument("-t", "--tsb", help="Sample name. Default parses from the file name", required =True)
    parser.add_argument("-b", "--build", help="Reference genome build [Default: hg38]", default="hg19")
    parser.add_argument("-p", "--protocol", help="Protocol used to generate annovar annotations [Default: wgEncodeGencodeBasicV19]", default="wgEncodeGencodeBasicV19")
    args = parser.parse_args()
    # print(args)


    meta = '\t'.join([str(x) for x in [args.tsb, args.build, "NA"]])

    with open(args.output, 'w') as handle:
        for output_line in read_annovar_file(args.input, meta, args.protocol):
            handle.write(output_line+"\n")
