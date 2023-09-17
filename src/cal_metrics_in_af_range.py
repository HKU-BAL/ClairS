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

from argparse import ArgumentParser
from collections import defaultdict

from shared.vcf import VcfReader

def cal_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return round(precision, 6), round(recall, 6), round(f1_score, 6)

def read_vcf(input_vcf_fn, ctg_name=None, filter_tag=None, keep_af=True, min_qual=None, max_qual = None):
    input_vcf_reader = VcfReader(vcf_fn=input_vcf_fn, ctg_name=ctg_name, show_ref=False, keep_row_str=True,
                                 skip_genotype=True, filter_tag=filter_tag, keep_af=keep_af, min_qual=min_qual,
                                 max_qual=max_qual)
    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict
    return input_variant_dict

def cal_metrics_in_af_range(args):
    low_af_path = os.path.join(args.low_af_path)
    af_range = args.af_range
    compare_vcf_output_dir = args.compare_vcf_output_dir

    if not os.path.exists(low_af_path):
        sys.exit('[ERROR] File not found: {}'.format(low_af_path))

    af_dict = defaultdict(float)
    af_range_list = af_range.split(',')

    output = open(low_af_path).read()
    for row in output.rstrip().split('\n'):
        if row == "":
            continue
        columns = row.rstrip().split(" ")
        ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt = columns[:6]
        tumor_af = float(tumor_alt) / (float(tumor_cov) if float(tumor_cov) > 0 else 1)
        af_dict[(ctg, int(pos))] = tumor_af

    fp_path = os.path.join(compare_vcf_output_dir, 'fp.vcf')
    tp_path = os.path.join(compare_vcf_output_dir, 'tp.vcf')
    fn_path = os.path.join(compare_vcf_output_dir, 'fn.vcf')
    for path in (fp_path, tp_path, fn_path):
        if not os.path.exists(path):
            sys.exit('[ERROR] File not found: {}'.format(path))

    fp = read_vcf(fp_path, keep_af=True, min_qual=args.min_qual)
    tp = read_vcf(tp_path, keep_af=True, min_qual=None)
    fn = read_vcf(fn_path, keep_af=False, min_qual=None)

    print("AF_range pre rec F1-score FP TP FN")

    #make sure include the variant with 1.0 AF
    total_fp, total_fn, total_tp = 0, 0, 0
    for af_idx, af in enumerate(af_range_list[:-1]):
        af = float(af)
        next_af = float(af_range_list[af_idx + 1])
        # include the max AF
        if next_af == 1.0:
            next_af = 1.01
        fp_in_af_key = [key for key in fp if float(fp[key].af) >= af and float(fp[key].af) < next_af]
        fp_in_af = len(fp_in_af_key)
        tp_in_af_key = [key for key in tp if key in af_dict and af_dict[key] >= af and af_dict[key] < next_af]
        tp_in_af = len(tp_in_af_key)
        if len(set(tp_in_af_key)) != len(tp_in_af_key):
            print("Not equal!!")
        # tp_dict[af] = set(tp_in_af_key)
        if args.min_qual is not None:

            tp_low_qual_in_af = sum([1 for key in tp if
                                     key in af_dict and af_dict[key] >= af and af_dict[key] < next_af and float(
                                         tp[key].qual) < args.min_qual])
        else:
            tp_low_qual_in_af = 0

        fn_in_af_key = [key for key in fn if key in af_dict and af_dict[key] >= af and af_dict[key] < next_af]
        fn_in_af = len(fn_in_af_key)

        pre, rec, f1 = cal_metrics(tp=tp_in_af - tp_low_qual_in_af, fp=fp_in_af, fn=fn_in_af + tp_low_qual_in_af)
        total_fp += fp_in_af
        total_tp += tp_in_af - tp_low_qual_in_af
        total_fn += fn_in_af + tp_low_qual_in_af
        if next_af == 1.01:
            next_af = 1.0
        output_str = "{}-{} {} {} {} {} {} {}".format(af, next_af, pre, rec, f1, fp_in_af, tp_in_af, fn_in_af)
        print(output_str)


def main():
    parser = ArgumentParser(description="Calculate metrics for various AF distribution")

    parser.add_argument('--low_af_path', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--af_range', type=str, default=None,
                        help="AF range, split by ','")

    parser.add_argument('--compare_vcf_output_dir', type=str, default=None,
                        help="The output dir of compare_vcf submodule")

    parser.add_argument('--min_qual', type=float, default=None,
                        help="")

    global args
    args = parser.parse_args()

    cal_metrics_in_af_range(args)


if __name__ == "__main__":
    main()
