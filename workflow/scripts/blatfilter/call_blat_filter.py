from blat_filter import RealignmentFilter
import argparse
import os

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser(description="Filter mutations with Blat.")
    parser.add_argument('--blat_executable', help='path to blat executable (jar file)', required=True)
    parser.add_argument('--bam', help='bam file', required=True)
    parser.add_argument('--maf', help='maf file. Lists mutations to filter.', required=True)
    parser.add_argument('--database', help='Reference 2bit file to blat sequences to.', required=True)
    parser.add_argument('--ooc', help='ooc specific to reference. Run blat -makeOoc', required=True)
    parser.add_argument(
        '--insert_size_metrics',
        help='Output metrics file from CollectMultipleMetrics or picard CollectInsertSizeMetrics',
        required=False
    )
    parser.add_argument(
        '--paired_ended',
        action="store_true"
    )
    parser.add_argument('--out_name', help='Output prefix for output maf files', default=None)
    parser.add_argument('--mapping_quality_threshold', type=int, default=5)
    parser.add_argument('--base_quality_threshold', type=int, default=5)
    parser.add_argument('--diff_score_concordant_threshold', type=int, default=3)
    parser.add_argument('--diff_score_discordant_threshold', type=int, default=-10)
    parser.add_argument('--concordant_fraction_threshold', type=float, default=0.75)
    parser.add_argument('--set_maf_col_manually', action='store_true')
    parser.set_defaults(set_maf_col_manually=False)

    parser.add_argument('--maf_chr_col', default='Chromosome')
    parser.add_argument('--maf_start_pos_col', default='Start_position')
    parser.add_argument('--maf_end_pos_col', default='End_position')
    parser.add_argument('--maf_ref_allele_col', default='Reference_Allele')
    parser.add_argument('--maf_alt_allele_col', default='Tumor_Seq_Allele2')
    parser.add_argument('--maf_t_ref_count_col', default='t_ref_count')
    parser.add_argument('--maf_t_alt_count_col', default='t_alt_count')
    parser.add_argument('--maf_chr_col_to_score_reads', default=None)
    parser.add_argument('--maf_start_pos_col_to_score_reads', default=None)
    parser.add_argument('--maf_pos_col_to_score_reads_empty_value', default=-1)
    parser.add_argument('--mate_insert_size_alpha', default=0.05)
    parser.add_argument('--stepper', default=0.05)
    args = parser.parse_args()
    # args_dict = vars(args)
    # print(args)
    args_dict = {
        "blat_executable": snakemake.params.blat_binary,
        "bam": snakemake.input.bam,
        "maf": snakemake.input.maf,
        "database": snakemake.input.database,
        "ooc": snakemake.input.occ,
        "paired_ended": True,
        "out_name": snakemake.params.output_prefix,
        "insert_size_metrics": snakemake.input.insert_size,
        "mapping_quality_threshold": 5,
        "base_quality_threshold": 5,
        "diff_score_concordant_threshold":3,
        "diff_score_discordant_threshold":10,
        "concordant_fraction_threshold": 0.75,
        "set_maf_col_manually": True,
        "maf_chr_col": "Chromosome",
        "maf_ref_allele_col": 'Reference_Allele',
        "maf_alt_allele_col": "Tumor_Seq_Allele2",
        "maf_chr_col_to_score_reads": "Chromosome",
        "maf_start_pos_col": "Start_Position",
        "maf_start_pos_col_to_score_reads": "Start_Position",
        "maf_end_pos_col": "End_Position",
        "maf_t_ref_count_col": "t_ref_count",
        "maf_t_alt_count_col": "t_alt_count",
        "maf_pos_col_to_score_reads_empty_value": -1,
        "mate_insert_size_alpha": 0.05,
        "genome": snakemake.input.genome,
    }
    stepper=snakemake.params.stepper
    print(stepper)
    print(args_dict['genome'])
    print(args_dict)
    os.makedirs(os.path.split(args_dict["out_name"])[0], exist_ok=True)
    RF = RealignmentFilter(**args_dict)

    # FIRST PASS

    print("\n**** first pass collecting blat queries ****\n")
    RF.collect_blat_queries(stepper=stepper)
    if RF._supporting_reads_df.empty:
        print("\n**** All varaints coming from anmolus reads rejecting all reads ****\n")
        RF.write_new_mafs_all_reject()
    else:

        print("\n**** first pass run blat aligner ****\n")
        RF.run_blat_aligner()

        print("\n**** first pass score supporting reads ****\n")
        RF.score_supporting_reads()

        print("\n**** first pass classify reads ****\n")
        RF.first_pass_process_read_scores()

        print("\n**** first pass filter variants ****\n")
        RF.first_pass_filter_judgements()

        # SECOND PASS
        if RF.paired_ended:
            print("\n**** second pass collecting blat mate queries ****\n")
            blat_mate_query_count = RF.collect_mate_blat_queries(stepper=stepper)

            if blat_mate_query_count > 0:

                print("\n**** second pass run blat aligner on mates ****\n")
                RF.run_blat_aligner(mate=True)

                print("\n**** second pass score mates ****\n")
                RF.score_supporting_reads(mate=True)

                print("\n**** second pass classify mates ****\n")
                RF.second_pass_process_read_scores()

                print("\n**** second pass filter variants ****\n")
                RF.second_pass_filter_judgements()

            else:
                print("\n**** second pass: no mate queries ****")
                print("**** second pass reject UNCERTAIN variants ****\n")
                RF.second_pass_reject_uncertain_variants()

        # WRITE RESULTS TO MAF FILES

        print("\n**** write new mafs ****\n")
        RF.write_new_mafs()
