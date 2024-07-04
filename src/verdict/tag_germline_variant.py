from argparse import ArgumentParser

import sys
import os

from argparse import ArgumentParser
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
from shared.vcf import VcfReader
import subprocess
from src.sort_vcf import compress_index_vcf
from numpy import *
from scipy.stats import binomtest

import random
from time import time

seed = int(time())
random.seed(seed)

def tag_germline_variant(args):

    input_vcf_fn = args.input_vcf_fn

    tumor_purity_path = args.tumor_purity_ploidy_output_file

    segment_path = args.tumor_cna_output_file
    if not os.path.exists(tumor_purity_path) or not os.path.exists(segment_path):
        return

    tumor_purity = float(open(tumor_purity_path).read().rstrip().split('\n')[1].split('\t')[1])
    input_vcf_reader = VcfReader(
            vcf_fn=input_vcf_fn,
            show_ref=True,
            keep_row_str=True,
            keep_af=True,
            skip_genotype=True,
            save_header=True,
            filter_tag='PASS'
    )
    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict

    cna_file = open(segment_path, 'r')

    seg_chr_list = []
    seg_start_list = []
    seg_end_list = []
    cn_major_list = []
    cn_minor_list = []
    for idx, cna in enumerate(cna_file.readlines()):
        if idx == 0:
            continue
        cna_columns = cna.strip().split('\t')
        seg_chr = str(cna_columns[1].strip('"'))
        seg_start = int(cna_columns[2])
        seg_end = int(cna_columns[3])
        cn_major = int(cna_columns[4])
        cn_minor = int(cna_columns[5])
        seg_chr_list.append(seg_chr)
        seg_start_list.append(seg_start)
        seg_end_list.append(seg_end)
        cn_major_list.append(cn_major)
        cn_minor_list.append(cn_minor)

    cna_file.close()
    high_purity_tresh = 0.95
    ALPHA = 0.01

    not_pass_count = 0
    for k, v in input_variant_dict.items():
        frequency = v.af
        columns = v.row_str.strip().split('\t')
        depth = int(columns[9].strip().split(':')[2])
        ctg_name = str(columns[0])
        pos = int(columns[1])
        p = M = C = nan
        segFound = 0
        clonality = 'NA'
        in_tumor = 'NA'
        SG_status = 'somatic'
        for i in range(len(seg_chr_list)):
            seg_chr = seg_chr_list[i]
            seg_start = seg_start_list[i]
            seg_end = seg_end_list[i]
            cn_major = cn_major_list[i]
            cn_minor = cn_minor_list[i]
            if seg_chr == ctg_name and seg_start <= pos <= seg_end:
                segFound = 1
                p = tumor_purity
                M = cn_minor
                C = cn_major + cn_minor
                if M == 0:
                    M = C - M
                P_G1 = P_G2 = P_S1 = P_S2 = 0
                AF_G1 = (p * M + 1 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                AF_S1 = (p * M + 0 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                P_G1 = binomtest(round(depth * frequency), depth, AF_G1).pvalue
                P_S1 = binomtest(round(depth * frequency), depth, AF_S1).pvalue
                if M != C - M:
                    AF_G2 = (p * (C - M) + 1 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                    P_G2 = binomtest(round(depth * frequency), depth, AF_G2).pvalue
                    if C - M != 0:
                        AF_S2 = (p * (C - M) + 0 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                        P_S2 = binomtest(round(depth * frequency), depth, AF_S2).pvalue
                    else:
                        AF_S2 = P_S2 = nan
                else:
                    AF_G2 = AF_S2 = P_G2 = P_S2 = nan

                max_prob_germline = max(P_G1, P_G2)
                max_prob_somatic = max(P_S1, P_S2)

                if max_prob_somatic == 0:
                    logodds = inf
                elif max_prob_germline == 0:
                    logodds = -inf
                else:
                    logodds = log10(max_prob_germline) - log10(max_prob_somatic)

                if frequency < 0.05 and 0.2 < p < 0.9:
                    SG_status = 'subclonal somatic'
                    columns[6] += ';SubclonalSomatic'

                elif isnan(C) or isnan(M):
                    SG_status = 'ambiguous_CNA_model'
                    # columns[6] += ';AmbiguousCNA'

                elif p > 0.95:
                    SG_status = 'nocall_purity>95%'
                    # columns[6] += ';NoCallPurity'

                elif frequency > 0.95:
                    SG_status = 'germline'
                    columns[6] += 'Germline'

                elif max_prob_germline > ALPHA and max_prob_somatic < ALPHA:
                    if logodds < 2:
                        SG_status = 'probable germline'
                        columns[6] += ';ProbableGermline'
                    else:
                        if frequency > 0.25:
                            SG_status = 'germline'
                            columns[6] += ';Germline'
                        else:
                            SG_status = 'probable germline'
                            columns[6] += ';ProbableGermline'

                elif max_prob_germline < ALPHA and max_prob_somatic > ALPHA:
                    if logodds > -2:
                        SG_status = 'probable somatic'
                        columns[6] += ';ProbableSomatic'
                    else:
                        SG_status = 'somatic'

                    if nanargmax([P_S1, P_S2]) == 1:
                        M = C - M

                elif max_prob_germline > ALPHA and max_prob_somatic > ALPHA:
                    SG_status = 'ambiguous_both_G_and_S'

                elif max_prob_germline < ALPHA and max_prob_somatic < ALPHA:

                    min_soma_EAF = min(AF_S1, AF_S2)
                    min_germ_EAF = min(AF_G1, AF_G2)

                    if p >= 0.3 and frequency < 0.25 and frequency < min_soma_EAF / 1.5 and min_soma_EAF <= min_germ_EAF:
                        SG_status = 'subclonal somatic'
                        columns[6] += ';SubclonalSomatic'

                    elif p >= 0.3 and frequency < 0.25 and frequency < min_germ_EAF / 2.0 and min_germ_EAF < min_soma_EAF:
                        SG_status = 'subclonal somatic'
                        columns[6] += ';SubclonalSomatic'

                    elif logodds < -5 and max_prob_somatic > 1e-10:
                        SG_status = 'somatic'
                        columns[6] += ';Somatic'

                        if nanargmax([P_S1, P_S2]) == 1:
                            M = C - M

                    elif logodds > 5 and max_prob_germline > 1e-4:
                        SG_status = 'germline'
                        columns[6] += 'Germline'

                    else:
                        SG_status = 'ambiguous_neither_G_nor_S'

                else:
                    SG_status = 'unknown'

                break


        v.row_str = '\t'.join(columns) + '\n'

    with open(args.output_fn, 'w') as f:
        #write header
        header = input_vcf_reader.header
        f.write(header)
        for k, v in input_variant_dict.items():
            f.write(v.row_str)

    if os.path.exists(args.output_fn):
        compress_index_vcf(args.output_fn)


def main():
    parser = ArgumentParser(description="")

    parser.add_argument('--tumor_purity_ploidy_output_file', type=str,
                        default=None,
                        help="Prefix of tumor allele count file")

    parser.add_argument('--tumor_cna_output_file', type=str,
                        default=None,
                        help="Prefix of normal allele count file")

    parser.add_argument('--input_vcf_fn', type=str,
                        default=None,
                        help="Prefix of 1kG alleles file")

    parser.add_argument('--output_fn', type=str,
                        default=None,
                        help="Output file path")

    global args
    args = parser.parse_args()
    tag_germline_variant(args)


if __name__ == "__main__":
    main()