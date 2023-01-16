# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import os
import logging
import shlex

from time import time
from argparse import ArgumentParser, SUPPRESS
from subprocess import run
from math import log, e
from collections import namedtuple

from shared.vcf import VcfWriter
from shared.utils import str2bool, subprocess_popen
import shared.param as param

logging.basicConfig(format='%(message)s', level=logging.INFO)

ACGT = 'ACGT'
Phred_Trans = (-10 * log(e, 10))

OutputConfig = namedtuple('OutputConfig', [
    'is_show_reference',
    'is_show_germline',
    'is_output_for_ensemble',
    'quality_score_for_pass',
    'tensor_fn',
    'input_probabilities',
    'pileup'
])
OutputUtilities = namedtuple('OutputUtilities', [
    'print_debug_message',
    'output',
    'output_header',
    'close_opened_files',
    'gen_output_file'
])


def filtration_value_from(quality_score_for_pass, quality_score, is_reference=False, is_germline=False):
    if is_reference:
        return 'RefCall'
    elif is_germline:
        return 'Germline'
    if quality_score_for_pass is None:
        return "PASS"
    if quality_score >= quality_score_for_pass:
        return "PASS"

    return "LowQual"


def quality_score_from(probability):
    return probability


def argmax(l):
    return max(zip(l, range(len(l))))[1]


def decode_acgt_count(alt_dict, ref_base=None, tumor_coverage=None):
    acgt_count = [0, 0, 0, 0]
    for idx, base in enumerate('ACGT'):
        acgt_count[idx] = alt_dict[base] if base in alt_dict else 0

    if ref_base is not None and tumor_coverage is not None:
        ref_base = ref_base[0].upper()
        if ref_base in 'ACGT':
            #update ref base count
            # acgt_count[]
            ref_idx = 'ACGT'.index(ref_base)
            acgt_count[ref_idx] = tumor_coverage - sum(acgt_count)

    AU, CU, GU, TU = acgt_count
    return AU, CU, GU, TU


def output_vcf_from_probability(
        chromosome,
        position,
        reference_base,
        normal_alt_info,
        tumor_alt_info,
        probabilities,
        output_config=None,
        vcf_writer=None,
):
    def decode_alt_info(alt_info):
        alt_info = alt_info.rstrip().split('-')
        read_depth = int(alt_info[0])  # alt_info
        indel_str = alt_info[1] if len(alt_info) > 1 else ''
        seqs = indel_str.split(' ')
        alt_info_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
        return alt_info_dict, read_depth

    normal_alt_info_dict, normal_read_depth = decode_alt_info(normal_alt_info)
    tumor_alt_info_dict, tumor_read_depth = decode_alt_info(tumor_alt_info)

    somatic_arg_index = param.somatic_arg_index
    alternate_base = reference_base
    arg_index = argmax(probabilities)
    is_reference = arg_index == 0
    is_germline = arg_index == 1
    is_tumor = arg_index == somatic_arg_index
    maximum_probability = probabilities[arg_index]

    def rank_somatic_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth):
        support_alt_dict = {}
        for tumor_alt, tumor_count in tumor_alt_info_dict.items():
            tumor_af = tumor_count / float(tumor_read_depth)
            normal_count = normal_alt_info_dict[tumor_alt] if tumor_alt in normal_alt_info_dict else 0
            normal_af = normal_count / float(normal_read_depth) if normal_read_depth > 0 else 0
            if tumor_af - normal_af > 0:
                support_alt_dict[tumor_alt] = tumor_af - normal_af
        if len(support_alt_dict) == 0:
            return "", 0, 0
        alt_type_list = sorted(support_alt_dict.items(), key=lambda x: x[1], reverse=True)
        best_match_alt = alt_type_list[0][0]
        tumor_supported_reads_count = tumor_alt_info_dict[best_match_alt]
        normal_supported_reads_count = normal_alt_info_dict[
            best_match_alt] if best_match_alt in normal_alt_info_dict else 0
        return best_match_alt, tumor_supported_reads_count, normal_supported_reads_count

    def rank_germline_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth):
        support_alt_dict = {}
        for tumor_alt, tumor_count in tumor_alt_info_dict.items():
            tumor_af = tumor_count / float(tumor_read_depth)
            normal_count = normal_alt_info_dict[tumor_alt] if tumor_alt in normal_alt_info_dict else 0
            normal_af = normal_count / float(normal_read_depth) if normal_read_depth > 0 else 0
            support_alt_dict[tumor_alt] = tumor_af + normal_af
        if len(support_alt_dict) == 0:
            return "", 0, 0
        alt_type_list = sorted(support_alt_dict.items(), key=lambda x: x[1], reverse=True)
        best_match_alt = alt_type_list[0][0]
        tumor_supported_reads_count = tumor_alt_info_dict[best_match_alt]
        normal_supported_reads_count = normal_alt_info_dict[
            best_match_alt] if best_match_alt in normal_alt_info_dict else 0
        return best_match_alt, tumor_supported_reads_count, normal_supported_reads_count

    if is_tumor:
        if tumor_read_depth <= 0:
            print("low tumor coverage")
            return
        best_match_alt, tumor_supported_reads_count, normal_supported_reads_count = rank_somatic_alt(
            tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth)

        if best_match_alt == "":
            return
        if best_match_alt[0] == 'X':
            alternate_base = best_match_alt[1]
            is_SNP = True
        elif best_match_alt[0] == 'I':
            alternate_base = best_match_alt[1:]
            is_INS = True
        elif best_match_alt[0] == 'D':
            alternate_base = reference_base
            reference_base += best_match_alt[1:]

    if is_germline and output_config.is_show_germline:
        best_match_alt, tumor_supported_reads_count, normal_supported_reads_count = rank_germline_alt(
            tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth)
        if best_match_alt == "":
            return
        if best_match_alt[0] == 'X':
            alternate_base = best_match_alt[1]
            is_SNP = True
        elif best_match_alt[0] == 'I':
            alternate_base = best_match_alt[1:]
            is_INS = True
        elif best_match_alt[0] == 'D':
            alternate_base = reference_base
            reference_base += best_match_alt[1:]

    if (not output_config.is_show_reference and is_reference) or (
            not is_reference and reference_base == alternate_base):
        return

    if (not output_config.is_show_germline and is_germline):
        return

    if reference_base is None or alternate_base is None:
        return

    # discard Indel
    if len(reference_base) > 1 or len(alternate_base) > 1:
        return

    def decode_alt_info(alt_info_dict, read_depth):
        alt_type_list = [{}, {}, {}]  # SNP I D
        snp_num, ins_num, del_num = 0, 0, 0
        for alt_type, count in alt_info_dict.items():
            count = int(count)
            if alt_type[0] == 'X':
                alt_type_list[0][alt_type[1]] = count
                snp_num += count
            elif alt_type[0] == 'I':
                alt_type_list[1][alt_type[1:]] = count
                ins_num += count
            elif alt_type[0] == 'D':
                alt_type_list[2][alt_type[1:]] = count
                del_num += count
        ref_num = max(0, read_depth - snp_num - ins_num - del_num)
        return alt_type_list, ref_num, snp_num, ins_num, del_num

    normal_alt_type_list, normal_ref_num, normal_snp_num, normal_ins_num, normal_del_num = decode_alt_info(
        alt_info_dict=normal_alt_info_dict, read_depth=normal_read_depth)
    tumor_alt_type_list, tumor_ref_num, tumor_snp_num, tumor_ins_num, tumor_del_num = decode_alt_info(
        alt_info_dict=tumor_alt_info_dict, read_depth=tumor_read_depth)

    if is_reference:
        normal_supported_reads_count = normal_ref_num
        tumor_supported_reads_count = tumor_ref_num
        alternate_base = "."

    normal_allele_frequency = min((normal_supported_reads_count / normal_read_depth) if normal_read_depth != 0 else 0.0,
                                  1.0)
    tumor_allele_frequency = min((tumor_supported_reads_count / tumor_read_depth) if tumor_read_depth != 0 else 0.0,
                                 1.0)

    # genotype string
    if is_reference:
        genotype_string = '0/0'
    elif is_germline:
        genotype_string = "0/1" if tumor_allele_frequency <= 0.5 else "0/1"
    else:
        genotype_string = "0/1" if tumor_allele_frequency < 1.0 else '1/1'
    # quality score
    quality_score = quality_score_from(maximum_probability)

    # filtration value
    filtration_value = filtration_value_from(
        quality_score_for_pass=output_config.quality_score_for_pass,
        quality_score=quality_score,
        is_reference=is_reference,
        is_germline=is_germline
    )

    information_string = "."

    AU, CU, GU, TU = decode_acgt_count(tumor_alt_type_list[0], reference_base, tumor_read_depth)

    vcf_writer.write_row(CHROM=chromosome,
                         POS=position,
                         REF=reference_base,
                         ALT=alternate_base,
                         QUAL=quality_score,
                         FILTER=filtration_value,
                         INFO=information_string,
                         GT=genotype_string,
                         DP=tumor_read_depth,
                         NDP=normal_read_depth,
                         AF=tumor_allele_frequency,
                         NAF=normal_allele_frequency,
                         AU=AU,
                         CU=CU,
                         GU=GU,
                         TU=TU)


def call_variants_from_probability(args):
    output_config = OutputConfig(
        is_show_reference=args.show_ref,
        is_show_germline=args.show_germline,
        is_output_for_ensemble=args.output_for_ensemble,
        quality_score_for_pass=args.qual,
        tensor_fn=args.tensor_fn,
        input_probabilities=args.input_probabilities,
        pileup=args.pileup
    )

    call_fn = args.call_fn
    if call_fn != "PIPE":
        call_dir = os.path.dirname(call_fn)
        if not os.path.exists(call_dir):
            output = run("mkdir -p {}".format(call_dir), shell=True)
    vcf_writer = VcfWriter(vcf_fn=args.call_fn,
                           ref_fn=args.ref_fn,
                           ctg_name=args.ctg_name,
                           show_ref_calls=args.show_ref,
                           sample_name=args.sample_name,
                           )

    logging.info("Calling somatic variants ...")
    variant_call_start_time = time()
    prediction_path = args.predict_fn

    if prediction_path != "PIPE":
        if not os.path.exists(prediction_path):
            print("[ERROR] Prediction path not found!")
            return
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, prediction_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    for row_id, row in enumerate(fo):
        row = row.rstrip().split('\t')
        chromosome, position, reference_base, normal_alt_info, tumor_alt_info, prediction = row[:6]
        probabilities = [float(item) for item in prediction.split()]
        output_vcf_from_probability(
            chromosome,
            position,
            reference_base,
            normal_alt_info,
            tumor_alt_info,
            probabilities,
            output_config=output_config,
            vcf_writer=vcf_writer)

    logging.info("Total time elapsed: %.2f s" % (time() - variant_call_start_time))

    vcf_writer.close()
    # remove file if on variant in output
    if os.path.exists(args.call_fn):
        vcf_file = open(args.call_fn, 'r').readlines()
        if not len(vcf_file):
            os.remove(args.call_fn)
        for row in vcf_file:
            if row[0] != '#':
                return
        logging.info("[INFO] No vcf output in file {}, remove.".format(args.call_fn))
        os.remove(args.call_fn)


def main():
    parser = ArgumentParser(description="Call variants using a trained model and tensors of candidate variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required if --gvcf is enabled")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=0,
                        help="If set, variants with >=QUAL will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to samtools, samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    parser.add_argument('--input_probabilities', action='store_true',
                        help="DEBUG: Use network probability outputs as input and generate variants from them")

    parser.add_argument('--output_probabilities', action='store_true',
                        help="DEBUG: Output the network probabilities of gt21, genotype, indel_length_1 and indel_length_2")

    # options for internal process control
    ## Use GPU for calling
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help=SUPPRESS)

    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Generating outputs for ensemble model calling
    parser.add_argument('--output_for_ensemble', action='store_true',
                        help=SUPPRESS)

    ## Use bin file from pytables to speed up calling.
    parser.add_argument('--is_from_tables', action='store_true',
                        help=SUPPRESS)

    ## Output reference calls
    parser.add_argument('--show_ref', action='store_true',
                        help=SUPPRESS)

    ## Output germline calls
    parser.add_argument('--show_germline', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    call_variants_from_probability(args)


if __name__ == "__main__":
    main()
