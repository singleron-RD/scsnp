#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd
import pysam
from sccore import utils


def extract_fields(record, samples):
    """
    extract fields from a snpeff annotated VCF record; ANN columns are 0-based
    CHROM, POS, REF, ALT, ANN[*].GENE(col 3), ANN[*].HGVS_C(col 9), ANN[*].HGVS_P(col 10)
    https://pcingola.github.io/SnpEff/snpsift/extractfields/
    """
    ann = record.info["ANN"][0].split("|")
    gene = ann[3]
    hgvs_c = ann[9]
    hgvs_p = ann[10]
    id = hgvs_p if hgvs_p else hgvs_c
    name = f"{gene}-{id}"

    genotypes = []
    for sample in samples:
        genotype = record.samples[sample]["GT"]
        if genotype == (None,):
            genotype_str = "NA"
        else:
            g1, g2 = genotype
            genotype_str = "/".join([str(g1), str(g2)])
        genotypes.append(genotype_str)

    return name, genotypes


def parse_vcf(vcf_fn):
    """
    Returns
    """
    vcf = pysam.VariantFile(vcf_fn)
    samples = list(vcf.header.samples)
    variant_dict = {}
    stats_dict = defaultdict(dict)
    for record in vcf:
        name, genotypes = extract_fields(record, samples)
        variant_dict[name] = genotypes
        for gt in ["0/0", "0/1", "1/1", "NA"]:
            stats_dict[name][gt] = genotypes.count(gt)
    df = pd.DataFrame(variant_dict, index=samples)
    return df, stats_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--vcf", required=True)
    args = parser.parse_args()

    gt_csv_fn = f"{args.sample}.GT.csv"
    stat_json_fn = f"{args.sample}.stats.json"

    df, stats_dict = parse_vcf(args.vcf)
    df.to_csv(gt_csv_fn)
    utils.write_json(stats_dict, stat_json_fn)


if __name__ == "__main__":
    main()
