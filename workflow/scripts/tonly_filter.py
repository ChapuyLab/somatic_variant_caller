import pandas as pd
import numpy as np
from scipy.stats import beta
import math
import sys
import argparse

# ====================== MAIN FUNCTION ======================
def purecn_predict_somatic(cns_file, vcf_file, out_vcf, purity, ploidy, contamination=0.0,
                           contamination_resolution=0.002, logodds_threshold=2.69,
                           tumor_sample="C-D-1145_tumor"):

    print(f"Running with purity={purity:.4f}, ploidy={ploidy:.2f}, contamination={contamination:.6f}")

    # 1. Load PureCN segments
    print("Reading PureCN segmentation...")
    cna = pd.read_csv(cns_file, sep="\t")
    cna = cna.rename(columns={"chromosome": "Chromosome", "start": "Start.bp", "end": "End.bp"})
    cna["NA"] = cna["cn1"].fillna(1).astype(float)
    cna["NB"] = cna["cn2"].fillna(cna["cn"] - cna["cn1"].fillna(0)).astype(float)
    cna["tau"] = cna["cn"].astype(float)
    cna["CCF_hat"] = 1.0
    cna["NCN"] = 2
    cna.loc[cna["Chromosome"] == "X", "NCN"] = 1
    cna.loc[cna["Chromosome"] == "Y", "NCN"] = 0

    # # 2. Parse only the variants we care about
    print("Parsing VCF (PASS / normal_artifact / germline only)...")
    maf = parse_vcf_to_maf_filtered(vcf_file, tumor_sample)
    print(f"\t{len(maf):,} variants kept for modeling")

    # 3. Run germline/somatic model
    print("Computing log-odds...")
    maf = germline_somatic_modeling(maf, cna, purity, contamination, ploidy,
                                    contamination_resolution, logodds_threshold)


    # 4. Write MAF outputs (unchanged)
    # out_base = f"{tumor_sample}.GermSomLogodds.purecn"
    # maf.to_csv(f"{out_base}.maf", sep="\t", index=False)
    # maf[maf["pass_tumor_only_filter"]].to_csv(f"{out_base}.pass.maf", sep="\t", index=False)

    # 5. NEW: Write annotated VCF with updated FILTER + VAF
    # out_vcf = f"{out_base}.vcf"
    print(f"Writing secondary VCF with updated FILTER and VAF → {out_vcf}")
    write_filtered_vcf(vcf_file, maf, out_vcf)

    print("Done!")
    # print(f"\tMAF full\t {out_base}.maf")
    # print(f"\tMAF PASS\t {out_base}.pass.maf")
    print(f"\tAnnotated VCF\t {out_vcf}")
    return maf


# ====================== VCF PARSER ======================
def parse_vcf_to_maf_filtered(vcf_file, tumor_sample):
    data = []
    kept = total = 0
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                    tumor_idx = header.index(tumor_sample)
                    format_idx = header.index("FORMAT")
                continue
            total += 1
            fields = line.strip().split("\t")
            filter_str = fields[6]
            ignore = [
                "blatReject", "clustered_events", "t_lod", "str_contraction",
                "read_position", "position", "fragment_length", "multiallelic", "clipping",
                "strand_artifact", "strand_bias", "slippage", "weak_evidence", "orientation",
                "haplotype", "normal_artifact", "panel_of_normals","map_qual", "contamination"
            ]

            if any(i in filter_str for i in ignore):
                continue
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            fmt_dict = dict(zip(fields[format_idx].split(":"), fields[tumor_idx].split(":")))
            ad = fmt_dict.get("AD", "0,0").split(",")
            t_ref = int(ad[0])
            t_alt = int(ad[1]) if len(ad) > 1 else 0
            data.append({
                "Chromosome": chrom,
                "Start_position": pos,
                "End_position": pos,
                "t_ref_count": t_ref,
                "t_alt_count": t_alt,
                "FILTER": filter_str
            })
            kept += 1
            if total % 1_000_000 == 0:
                print(f"   parsed {total:,} → {kept:,} kept")
    print(f"   VCF scan done: {kept:,} kept")
    print(pd.DataFrame(data)["FILTER"].unique())
    return pd.DataFrame(data)


# ====================== CORE MODEL (unchanged) ======================
def germline_somatic_modeling(mutation_data, CNA_data, purity, contamination,
                              ploidy, CONTAMINATION_RESOLUTION, logodds_threshold):
    contamination = math.sqrt(contamination**2 + CONTAMINATION_RESOLUTION**2)
    mutation_data = assign_scna_state(mutation_data, CNA_data)
    NALT = mutation_data["t_alt_count"].astype(float)
    NREF = mutation_data["t_ref_count"].astype(float)
    i_tumor_f = NALT / (NALT + NREF + 1e-9)
    SCNA_NA = mutation_data["SCNA_NA"]
    SCNA_NB = mutation_data["SCNA_NB"]
    SCNA_NCN = mutation_data["SCNA_NCN"]
    SCNA_CCF = mutation_data["SCNA_CCF_hat"]
    DNA = ((1 - purity) * SCNA_NCN +
           purity * (1 - SCNA_CCF) * SCNA_NCN +
           purity * SCNA_CCF * (SCNA_NA + SCNA_NB)) + \
          (contamination * (ploidy * purity + (1 - purity) * 2) / (1 - contamination))
    pNULL = beta.pdf(i_tumor_f, NALT + 1, NREF + 1)
    GA = (1 - purity) + purity * (1 - SCNA_CCF) + purity * SCNA_CCF * SCNA_NA
    GB = (1 - purity) + purity * (1 - SCNA_CCF) + purity * SCNA_CCF * SCNA_NB
    G_HOM = (1 - purity) * SCNA_NCN + purity * (1 - SCNA_CCF) * SCNA_NCN + purity * SCNA_CCF * (SCNA_NA + SCNA_NB)
    LgA = beta.pdf(GA / DNA, NALT + 1, NREF + 1)
    LgB = beta.pdf(GB / DNA, NALT + 1, NREF + 1)
    LgHom = beta.pdf(G_HOM / DNA, NALT + 1, NREF + 1)
    S1 = purity * SCNA_CCF * (1 - SCNA_CCF)
    S2 = purity * SCNA_CCF * SCNA_CCF
    S3 = purity * SCNA_CCF * SCNA_CCF * SCNA_NA
    S4 = purity * SCNA_CCF * SCNA_CCF * SCNA_NB
    S5 = purity * SCNA_CCF * (1 - SCNA_CCF) + purity * SCNA_CCF * SCNA_CCF * SCNA_NA
    S6 = purity * SCNA_CCF * (1 - SCNA_CCF) + purity * SCNA_CCF * SCNA_CCF * SCNA_NB
    S_L1 = beta.pdf(S1 / DNA, NALT + 1, NREF + 1)
    S_L2 = beta.pdf(S2 / DNA, NALT + 1, NREF + 1)
    S_L3 = beta.pdf(S3 / DNA, NALT + 1, NREF + 1)
    S_L4 = beta.pdf(S4 / DNA, NALT + 1, NREF + 1)
    S_L5 = beta.pdf(S5 / DNA, NALT + 1, NREF + 1)
    S_L6 = beta.pdf(S6 / DNA, NALT + 1, NREF + 1)
    mutation_data["logodds_Germline_Somatic"] = np.log10(
        np.maximum.reduce([LgA, LgB, LgHom]) /
        np.maximum.reduce([S_L1, S_L2, S_L3, S_L4, S_L5, S_L6])
    )
    gs_cutoff = logodds_threshold * purity
    mutation_data["GermSom_logodds_cutoff"] = np.maximum(0, gs_cutoff)
    mutation_data["pass_tumor_only_filter"] = (
        mutation_data["logodds_Germline_Somatic"] <= mutation_data["GermSom_logodds_cutoff"]
    )
    mutation_data.loc[mutation_data["logodds_Germline_Somatic"].isna(), "pass_tumor_only_filter"] = True
    return mutation_data


def assign_scna_state(mutation_data, CNA_data):
    mutation_data = mutation_data.copy()
    mutation_data[["SCNA_NA", "SCNA_NB", "SCNA_CCF_hat", "SCNA_NCN"]] = np.nan
    for _, row in CNA_data.iterrows():
        mask = ((mutation_data["Chromosome"] == str(row["Chromosome"])) &
                (mutation_data["Start_position"] >= row["Start.bp"]) &
                (mutation_data["End_position"] <= row["End.bp"]))
        if mask.any():
            mutation_data.loc[mask, ["SCNA_NA", "SCNA_NB", "SCNA_CCF_hat", "SCNA_NCN"]] = [
                row["NA"], row["NB"], row["CCF_hat"], row["NCN"]
            ]
    return mutation_data


# ====================== NEW: WRITE ANNOTATED VCF ======================
def write_filtered_vcf(original_vcf_path, maf, output_vcf_path):
    # Build fast lookup: (chrom, pos) → (is_somatic, original_filter, vaf)
    lookup = {}
    for _, row in maf.iterrows():
        key = (row["Chromosome"], row["Start_position"])
        vaf = row["t_alt_count"] / (row["t_alt_count"] + row["t_ref_count"] + 1e-9)
        lookup[key] = (bool(row["pass_tumor_only_filter"]), row["FILTER"], round(vaf, 6))

    with open(original_vcf_path) as fin, open(output_vcf_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    fout.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Tumor variant allele fraction (t_alt/(t_alt+t_ref))">\n')
                    fout.write('##FILTER=<ID=TonlyFilter,Description="Failed tumor-only germline/somatic logodds filter">\n')
                    fout.write(line)
                    continue
                fout.write(line)
                continue

            # Variant line
            fields = line.strip().split("\t")
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            key = (chrom, pos)

            if key in lookup:
                is_somatic, orig_filter, vaf = lookup[key]
                if orig_filter == "PASS":                    # only touch original PASS variants
                    fields[6] = "PASS" if is_somatic else "TonlyFilter"
                # else: keep original FILTER (normal_artifact, germline, etc.)

                # Add VAF to INFO
                info = fields[7]
                if info == ".":
                    info = f"VAF={vaf}"
                else:
                    info += f";VAF={vaf}"
                fields[7] = info

            fout.write("\t".join(fields) + "\n")


# ====================== CLI ======================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PureCN-like tumor-only filter + annotated VCF output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("cns_file", help="PureCN .cns file")
    parser.add_argument("vcf_file", help="Mutect2 VCF")
    parser.add_argument("outut_vcf_file", help="Output VCF")
    parser.add_argument("--purity", type=float, required=True, help="Tumor purity (0 < purity ≤ 1)")
    parser.add_argument("--ploidy", type=float, default=2.0, help="Tumor ploidy")
    parser.add_argument("--contamination", type=float, default=0.0, help="Contamination fraction")
    parser.add_argument("--contamination-table", type=str, default=None,
                        help="Mutect2 contamination.table – will override --contamination")
    parser.add_argument("--tumor-sample", default="C-D-1145_tumor", help="Tumor sample name in VCF")
    parser.add_argument("--logodds-threshold", type=float, default=2.69, help="Logodds cutoff coefficient")

    args = parser.parse_args()

    if args.contamination_table:
        cont_df = pd.read_csv(args.contamination_table, sep="\t")
        args.contamination = float(cont_df.iloc[0]["contamination"])
        print(f"\t Read contamination = {args.contamination:.6f} from table")

    if not 0 < args.purity <= 1:
        parser.error("Purity must be >0 and ≤1")

    purecn_predict_somatic(
        cns_file=args.cns_file,
        vcf_file=args.vcf_file,
        out_vcf=args.outut_vcf_file,
        purity=args.purity,
        ploidy=args.ploidy,
        contamination=args.contamination,
        logodds_threshold=args.logodds_threshold,
        tumor_sample=args.tumor_sample
    )
