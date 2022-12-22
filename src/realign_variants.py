import sys
import os
import subprocess
import concurrent.futures

from collections import Counter
from argparse import ArgumentParser
from collections import defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter
from shared.utils import str2bool

file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
main_entry = os.path.join(file_directory, "clair-somatic.py")


def get_base_list(columns):
    pileup_bases = columns[4]

    base_idx = 0
    base_list = []
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
    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list


def extract_base(POS):
    pos = POS.pos
    ref_base = POS.reference_bases
    alt_base = POS.alternate_bases[0]

    bam_fn = args.bam_fn
    ref_fn = args.ref_fn
    samtools = args.samtools
    ctg_name = args.ctg_name if args.ctg_name is not None else POS.ctg_name
    min_mq = args.min_mq
    min_bq = args.min_bq
    python = args.python
    qual = float(POS.qual) if POS.qual is not None else None
    if POS.extra_infos is False or (qual is not None and qual >= 0.95):
        return ctg_name, pos, True, (-1, -1, -1, -1)

    ctg_range = "{}:{}-{}".format(ctg_name, pos, pos)
    samtools_command = "{} mpileup {} --min-MQ {} --min-BQ {} --excl-flags 2316 -r {}".format(samtools,
                                                                                              bam_fn,
                                                                                              min_mq,
                                                                                              min_bq,
                                                                                              ctg_range)

    output = subprocess.run(samtools_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output = output.stdout.rstrip()

    columns = output.split('\t')
    if len(columns) < 4:
        return ctg_name, pos, True, (-1, -1, -1, -1)
    base_counter, base_list = get_base_list(columns)

    realign_command = "{} {} realign_reads --pos {} --ctg_name {} --bam_fn {} --ref_fn {} --samtools {}".format(python,
                                                                                                                main_entry,
                                                                                                                pos,
                                                                                                                ctg_name,
                                                                                                                bam_fn,
                                                                                                                ref_fn,
                                                                                                                samtools)

    samtools_mpileup_command = "{} mpileup - --reverse-del --min-MQ {} --min-BQ {} --excl-flags 2316 | grep -w {}".format(
                                                                                                                samtools,
                                                                                                                min_mq,
                                                                                                                min_bq,
                                                                                                                pos)
    realign_command += " | " + samtools_mpileup_command
    realign_output = subprocess.run(realign_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    realign_output = realign_output.stdout.rstrip()
    columns = realign_output.split('\t')
    if len(columns) < 4:
        return ctg_name, pos, True, (-1, -1, -1, -1)

    realign_base_counter, realign_base_list = get_base_list(columns)

    raw_depth = len(base_list)
    realign_depth = len(realign_base_list)
    raw_support_read_num = base_counter[alt_base]
    realign_support_read_num = realign_base_counter[alt_base]

    pass_realign_filter = True
    if raw_support_read_num / float(
            raw_depth) > realign_support_read_num / realign_depth and realign_support_read_num < raw_support_read_num:
        pass_realign_filter = False
    return ctg_name, pos, pass_realign_filter, (
    raw_support_read_num, raw_depth, realign_support_read_num, realign_depth)


def realign_variants(args):
    ctg_name = args.ctg_name
    threads = args.threads
    output_dir = args.output_dir
    enable_realignment = args.enable_realignment
    qual_cut_off = args.qual if args.qual is not None else param.qual_dict['ilmn']

    pileup_output_vcf_fn = os.path.join(output_dir, "pileup_filter.vcf")
    fa_output_vcf_fn = os.path.join(output_dir, "full_alignment_filter.vcf")
    if not enable_realignment:
        subprocess.run("ln -sf {} {}".format(args.pileup_vcf_fn, pileup_output_vcf_fn), shell=True)
        subprocess.run("ln -sf {} {}".format(args.full_alignment_vcf_fn, fa_output_vcf_fn), shell=True)
        return

    p_reader = VcfReader(vcf_fn=args.pileup_vcf_fn,
                         ctg_name=ctg_name,
                         show_ref=True,
                         keep_row_str=True,
                         filter_tag=None,
                         save_header=True)
    p_reader.read_vcf()
    p_input_variant_dict = p_reader.variant_dict

    fa_reader = VcfReader(vcf_fn=args.full_alignment_vcf_fn,
                          ctg_name=ctg_name,
                          show_ref=True,
                          keep_row_str=True,
                          filter_tag=None,
                          save_header=True)
    fa_reader.read_vcf()
    fa_input_variant_dict = fa_reader.variant_dict

    p_vcf_writer = VcfWriter(vcf_fn=pileup_output_vcf_fn,
                             ctg_name=ctg_name,
                             ref_fn=args.ref_fn,
                             show_ref_calls=True)

    f_vcf_writer = VcfWriter(vcf_fn=fa_output_vcf_fn,
                             ctg_name=ctg_name,
                             ref_fn=args.ref_fn,
                             show_ref_calls=True)

    for k, v in fa_input_variant_dict.items():
        if v.filter != "PASS":
            v.extra_infos = False
            continue
        if qual_cut_off is not None and float(v.qual) <= qual_cut_off:
            if k not in p_input_variant_dict:
                v.extra_infos = False

    total_num = 0
    realign_fail_pos_set = set()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as exec:
        for result in exec.map(extract_base, list(fa_input_variant_dict.values())):
            contig, pos, pass_realign_filter = result[:3]
            if pass_realign_filter is False:
                realign_fail_pos_set.add((contig, pos))
            total_num += 1
            if total_num > 0 and total_num % 1000 == 0:
                print("[INFO] Processing in {}, total processed positions: {}".format(contig, total_num))

    #write output
    for k, v in p_input_variant_dict.items():
        pos = k if ctg_name is not None else k[1]
        contig = ctg_name if ctg_name is not None else k[0]
        row_str = v.row_str.rstrip()
        if (contig, pos) in realign_fail_pos_set:
            row_str = row_str.replace("PASS", "LowQual")

        p_vcf_writer.vcf_writer.write(row_str + '\n')

    for k, v in fa_input_variant_dict.items():
        pos = k if ctg_name is not None else k[1]
        contig = ctg_name if ctg_name is not None else k[0]
        row_str = v.row_str.rstrip()
        if (contig, pos) in realign_fail_pos_set:
            row_str = row_str.replace("PASS", "LowQual")

        f_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()
    f_vcf_writer.close()

    if ctg_name is not None:
        print("[INFO] {}: {} called variant filtered by short-read realignment".format(ctg_name, len(realign_fail_pos_set)))
    else:
        print("[INFO] {} called variant filtered by realignment".format(len(realign_fail_pos_set)))


def main():
    parser = ArgumentParser(description="Reads realignment")

    parser.add_argument('--bam_fn', type=str, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF file input")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Full-alignment VCF file input")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output vcf directory")

    # options for advanced users
    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    # options for debug purpose
    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--threads', type=int, default=1,
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--pos', type=int, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")

    parser.add_argument('--enable_realignment', type=str2bool, default=True,
                        help="Enable realignment for illumina calling, default: enable")

    parser.add_argument('--qual', type=float, default=None,
                        help="Maximum QUAL to realign a variant")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    global args
    args = parser.parse_args()

    realign_variants(args)


if __name__ == "__main__":
    main()
