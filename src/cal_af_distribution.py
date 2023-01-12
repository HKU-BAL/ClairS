import sys
import os
import subprocess
import concurrent.futures

from collections import Counter
from argparse import ArgumentParser
from collections import defaultdict

from shared.vcf import VcfReader
import shared.param as param
from shared.utils import str2bool

def get_base_list(columns, args=None):
    if len(columns) < 5:
        return Counter(), []
    min_bq_cut = args.min_bq_cut if args is not None else 0
    pileup_bases = columns[4]
    base_idx = 0
    base_list = []
    bq_list = [ord(qual) - 33 for qual in columns[5]]
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
        # skip $, the end of read
        base_idx += 1
    upper_base_counter = Counter([''.join(item).upper() for item, bq in zip(base_list, bq_list) if bq >= min_bq_cut])
    return upper_base_counter, base_list


def extract_base(POS):

    pos = POS.pos
    ref_base = POS.reference_bases
    alt_base = POS.alternate_bases[0].upper()
    ctg_name = POS.ctg_name

    normal_samtools_command_with_region = normal_samtools_command + ' -r {}:{}-{}'.format(ctg_name, pos, pos)
    tumor_samtools_command_with_region = tumor_samtools_command + ' -r {}:{}-{}'.format(ctg_name, pos, pos)
    #normal
    output = subprocess.run(normal_samtools_command_with_region, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    output = output.stdout.rstrip()
    columns = output.split('\t')
    base_counter, base_list = get_base_list(columns)

    #tumor
    tumor_output = subprocess.run(tumor_samtools_command_with_region, shell=True, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, universal_newlines=True)
    tumor_output = tumor_output.stdout.rstrip()
    columns = tumor_output.split('\t')
    tumor_base_counter, tumor_base_list = get_base_list(columns)

    HAP_LIST = [0, 0, 0]
    ALL_HAP_LIST = [0, 0, 0]
    if len(columns) >= 7:
        phasing_info = columns[6].split(',')
        for hap_idx, (b, hap) in enumerate(zip(tumor_base_list, phasing_info)):
            if hap not in '12':
                hap = 0
            ALL_HAP_LIST[int(hap)] += 1
            if ''.join(b).upper() == alt_base:
                HAP_LIST[int(hap)] += 1

    HAP_LIST = " ".join([str(i) for i in HAP_LIST])
    ALL_HAP_LIST = " ".join([str(i) for i in ALL_HAP_LIST])
    return [ctg_name, pos, len(base_list), len(tumor_base_list), base_counter[alt_base], tumor_base_counter[alt_base],
            HAP_LIST, ALL_HAP_LIST]


class INFO():
    def __init__(self):
        self.normal_base_counter = None
        self.normal_depth = None
        self.tumor_base_counter = None
        self.tumor_depth = None
        self.alt_base = None


def parser_info(row):
    columns = row.split('\t')
    info_list = columns[8].split(':')
    FORMAT = columns[9].split(':')
    DP_index = info_list.index("DP")
    NDP_index = info_list.index('NDP')
    AF_index = info_list.index('AF')
    NAF_index = info_list.index('NAF')
    tumor_depth = int(FORMAT[DP_index])
    normal_depth = int(FORMAT[NDP_index])
    tumor_alt_depth = round(float(FORMAT[AF_index]) * tumor_depth)
    normal_alt_depth = round(float(FORMAT[NAF_index]) * normal_depth)

    return columns[0], columns[1], normal_depth, tumor_depth, normal_alt_depth, tumor_alt_depth


def cal_af(args, truth_variant_dict=None, input_variant_dict=None):
    ctg_name = args.ctg_name
    output_path = args.output_path

    if truth_variant_dict is None:
        truth_vcf_fn = args.truth_vcf_fn
        vcf_reader = VcfReader(vcf_fn=truth_vcf_fn,
                               ctg_name=ctg_name,
                               show_ref=False,
                               keep_row_str=True,
                               filter_tag=args.truth_filter_tag)
        vcf_reader.read_vcf()
        truth_variant_dict = vcf_reader.variant_dict

    if input_variant_dict is None:
        input_vcf_fn = args.input_vcf_fn
        vcf_reader = VcfReader(vcf_fn=input_vcf_fn,
                               ctg_name=ctg_name,
                               show_ref=False,
                               keep_row_str=True,
                               filter_tag=args.input_filter_tag)
        vcf_reader.read_vcf()
        input_variant_dict = vcf_reader.variant_dict

    results_dict = defaultdict()
    variant_dict = defaultdict()

    if output_path is not None:
        output_file = open(output_path, 'w')
    for k, v in truth_variant_dict.items():
        if k not in input_variant_dict:
            variant_dict[k] = v
        else:
            result = parser_info(input_variant_dict[k].row_str)
            if output_path is not None:
                output_file.write(' '.join(str(item) for item in result) + '\n')
            else:
                key = k if args.ctg_name is None else (args.ctg_name, k)
                results_dict[key] = result

    min_mq = param.min_mq
    min_bq = param.min_bq

    phasing_option = "--output-extra HP " if args.phase_output else " "
    samtools_command = "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 {} ".format(args.samtools,
                                                                                          min_mq,
                                                                                          min_bq,
                                                                                          phasing_option)

    global normal_samtools_command, tumor_samtools_command
    normal_samtools_command = samtools_command + args.normal_bam_fn
    tumor_samtools_command = samtools_command + args.tumor_bam_fn

    total_num = 0
    print("[INFO] Total truth need to calculate AF: {}".format(len(variant_dict)))

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as exec:
        for result in exec.map(extract_base, list(variant_dict.values())):
            total_num += 1
            if total_num % 1000 == 0 and total_num > 0:
                print("[INFO] Total processed positions: {}".format(total_num))
            if result is not None:
                ctg_name, pos, normal_depth, tumor_depth, normal_alt_depth, tumor_alt_depth, HAP_LIST, ALL_HAP_LIST = result
                k = (ctg_name, int(pos))
                results_dict[k] = ctg_name, pos, normal_depth, tumor_depth, normal_alt_depth, tumor_alt_depth, HAP_LIST, ALL_HAP_LIST
                if output_path is not None:
                    output_file.write(' '.join(str(item) for item in result) + '\n')

    if output_path is not None:
        output_file.close()
        return

    return results_dict


def main():
    parser = ArgumentParser(description="Calculate AF distribution in tumor and normal BAMs")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input, required")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input vcf filename")

    parser.add_argument('--truth_vcf_fn', type=str, default="PASS;HighConf,PASS;MedConf",
                        help="Truth vcf filename")

    parser.add_argument('--input_filter_tag', type=str, default=None,
                        help="Filter tag for the input VCF")

    parser.add_argument('--truth_filter_tag', type=str, default=None,
                        help="Filter tag for the truth VCF")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--output_path', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--min_bq_cut', type=int, default=0,
                        help="Output VCF filename, required")

    parser.add_argument('--phase_output', type=str2bool, default=False,
                        help="Output phasing INFO")

    parser.add_argument('--pos', type=int, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--threads', type=int, default=4,
                        help="Max #threads to be used")

    global args
    args = parser.parse_args()

    cal_af(args)


if __name__ == "__main__":
    main()
