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

import os
import sys
import subprocess
from argparse import ArgumentParser, SUPPRESS

from collections import defaultdict
from shared.utils import str2bool, str_none
from shared.vcf import VcfReader, VcfWriter
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import file_path_from
from src.cal_af_distribution import cal_af

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

major_contigs = {"chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]}.union(
    {str(a) for a in list(range(1, 23)) + ["X", "Y"]})


def cal_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return round(precision, 4), round(recall, 4), round(f1_score, 4)


def output_best_cut_off(fp_qual_dict, tp_qual_dict, fn_count, use_int_cut_off=True, add_tp_fn=False):
    results = []
    if use_int_cut_off:
        qual_list = set([int(q) for q in list(fp_qual_dict.values()) + list(tp_qual_dict.values())])
    else:
        qual_list = [item / 100.0 for item in range(0, 101)]

    for qual in qual_list:
        fp_snv = sum([1 for k, v in fp_qual_dict.items() if v >= qual])
        tp_snv = sum([1 for k, v in tp_qual_dict.items() if v >= qual])
        fn_snv = fn_count + len(tp_qual_dict) - tp_snv
        snv_pre, snv_rec, snv_f1 = cal_metrics(tp=tp_snv, fp=fp_snv, fn=fn_snv)
        tp_fn = tp_snv + fn_snv
        results.append([qual, snv_pre, snv_rec, snv_f1, tp_snv, fp_snv, fn_snv, tp_fn])

    results = sorted(results, key=lambda x: x[3], reverse=True)
    return results

def compare_vcf(args):
    """
    Follow how som.py works
    ## https://github.com/Illumina/hap.py/blob/master/doc/sompy.md
    """
    output_fn = args.output_fn
    output_dir = args.output_dir
    truth_vcf_fn = args.truth_vcf_fn
    input_vcf_fn = args.input_vcf_fn
    bed_fn = args.bed_fn
    high_confident_only = args.high_confident_only
    ctg_name = args.ctg_name
    skip_genotyping = args.skip_genotyping
    input_filter_tag = args.input_filter_tag
    truth_filter_tag = args.truth_filter_tag
    discard_fn_out_of_fp_bed = args.discard_fn_out_of_fp_bed
    benchmark_indel = args.benchmark_indel

    fp_bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=ctg_name)
    strat_bed_tree_list = []

    if args.strat_bed_fn is not None and ',' in args.strat_bed_fn:
        for strat_bed_fn in args.strat_bed_fn.split(','):
            strat_bed_tree_list.append(bed_tree_from(bed_file_path=strat_bed_fn, contig_name=ctg_name))
    elif args.strat_bed_fn is not None:
        strat_bed_tree_list = [bed_tree_from(bed_file_path=args.strat_bed_fn, contig_name=ctg_name)]

    truth_vcf_fn = file_path_from(file_name=truth_vcf_fn, exit_on_not_found=True, allow_none=False)
    input_vcf_fn = file_path_from(file_name=input_vcf_fn, exit_on_not_found=True, allow_none=False)

    truth_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn,
                                 ctg_name=ctg_name,
                                 ctg_start=args.ctg_start,
                                 ctg_end=args.ctg_end,
                                 show_ref=False,
                                 keep_row_str=True,
                                 skip_genotype=skip_genotyping,
                                 filter_tag=truth_filter_tag)
    truth_vcf_reader.read_vcf()
    truth_variant_dict = truth_vcf_reader.variant_dict

    input_vcf_reader = VcfReader(vcf_fn=input_vcf_fn,
                                 ctg_name=ctg_name,
                                 ctg_start=args.ctg_start,
                                 ctg_end=args.ctg_end,
                                 show_ref=False,
                                 keep_row_str=True,
                                 skip_genotype=skip_genotyping,
                                 filter_tag=input_filter_tag,
                                 keep_af=True,
                                 min_qual=args.min_qual,
                                 max_qual=args.max_qual,
                                 naf_filter=args.naf_filter,
                                 discard_indel=False if benchmark_indel else True)
    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict

    input_out_of_bed = 0
    truth_out_of_bed = 0
    low_bq_count = 0
    low_qual_truth = set()
    low_af_truth = set()
    input_out_of_strat_bed, truth_out_of_strat_bed = 0, 0

    discard_low_af = args.min_af is not None
    no_normal_count, no_tumor_alt_count, low_tumor_count, high_normal_af, low_af_count = 0, 0, 0, 0, 0
    if discard_low_af:
        if args.low_af_path is None or not os.path.exists(args.low_af_path):
            if args.normal_bam_fn is None or args.tumor_bam_fn is None or (not (os.path.exists(args.normal_bam_fn) and os.path.exists(args.tumor_bam_fn))):
                sys.exit("[ERROR] Pls input --normal_bam_fn and --tumor_bam_fn for calculation")
            result_dict = cal_af(args, truth_variant_dict, input_variant_dict)
        else:
            result_dict = defaultdict()
            fp = open(args.low_af_path).readlines()
            fp = [item.rstrip().split(' ') for item in fp]
            for row in fp:
                ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, *hap_info = row
                result_dict[ctg, int(pos)] = ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, hap_info

        for k, v in result_dict.items():
            ctg_name, pos, normal_cov, tumor_cov, normal_alt, tumor_alt = v[:6]
            key = int(pos) if args.ctg_name is not None else (ctg_name, int(pos))
            if int(normal_cov) <= args.min_coverage:
                low_af_truth.add(key)
                no_normal_count += 1
            elif int(tumor_alt) == 0 or int(tumor_cov) == 0:
                low_af_truth.add(key)
                no_tumor_alt_count += 1
            elif int(tumor_alt) / float(tumor_cov) <= args.min_af or int(tumor_alt) <= args.min_alt_coverage:
                low_tumor_count += 1
                low_af_truth.add(key)

        for k in list(input_variant_dict.keys()):
            if float(input_variant_dict[k].af) <= args.min_af:
                del input_variant_dict[k]
                low_af_count += 1

    phasable_count = 0
    non_phasable_count = 0
    if args.validate_phase_only is not None:
        unphase = 0
        phase_dict = defaultdict()
        if args.phase_output is not None:
            for item in open(args.phase_output).readlines():
                ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, *hap_info = item.rstrip().split(' ')
                phase_dict[ctg, int(pos)] = (ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, hap_info)
        else:
            phase_dict = result_dict
        for item in phase_dict.values():
            ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, hap_info = item
            hp0, hp1, hp2, all_hp0, all_hp1, all_hp2 = [int(i) for i in hap_info]
            key = int(pos) if args.ctg_name is not None else (ctg, int(pos))
            if key in input_variant_dict:
                continue

            phaseable = all_hp1 * all_hp2 > 0 and hp1 * hp2 == 0 and (
                        int(hp1) > args.min_alt_coverage or int(hp2) > args.min_alt_coverage)
            # phaseable 1 non phaseable 2
            if int(args.validate_phase_only) == 1 and not phaseable:
                low_af_truth.add(key)
                continue
            if int(args.validate_phase_only) == 2 and phaseable:
                low_af_truth.add(key)
                continue

        for k, v in input_variant_dict.items():
            columns = v.row_str.rstrip().split('\t')
            phaseable = columns[7] == 'H'
            if phaseable:
                phasable_count += 1
            else:
                non_phasable_count += 1

            if int(args.validate_phase_only) == 1 and not phaseable:
                low_af_truth.add(k)
                continue
            if int(args.validate_phase_only) == 2 and phaseable:
                low_af_truth.add(k)
                continue

        print("[INFO] Fail or non-phasable:", non_phasable_count, 'Low ALT base phase count :', unphase, "Phasable:",
              phasable_count)

    if high_confident_only:
        for key in list(truth_variant_dict.keys()):
            row = truth_variant_dict[key].row_str
            if "PASS;HighConf" not in row:
                low_qual_truth.add(key)

    if output_fn:
        output_file = open(output_fn, 'w')
    else:
        output_file = None

    for key in list(input_variant_dict.keys()):
        pos = key if args.ctg_name is not None else key[1]
        contig = args.ctg_name if args.ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                                contig_name=contig,
                                                                region_start=pos - 1,
                                                                region_end=pos)
        if not pass_bed_region:
            del input_variant_dict[key]
            input_out_of_bed += 1
            continue

        pass_straed_region = len(strat_bed_tree_list) == 0 or sum([1 if is_region_in(tree=strat_bed_tree,
                                                                                     contig_name=contig,
                                                                                     region_start=pos - 1,
                                                                                     region_end=pos) else 0
                                                                   for strat_bed_tree in strat_bed_tree_list]) == len(
            strat_bed_tree_list)

        if not pass_straed_region and key in input_variant_dict:
            del input_variant_dict[key]
            input_out_of_strat_bed += 1
            continue

        if high_confident_only and key in low_qual_truth:
            continue

        if benchmark_indel:
            ref_base, alt_base = input_variant_dict[key].reference_bases, input_variant_dict[key].alternate_bases[0]
            if len(ref_base) == 1 and len(alt_base) == 1:
                del input_variant_dict[key]

    for key in list(truth_variant_dict.keys()):
        pos = key if args.ctg_name is not None else key[1]
        contig = args.ctg_name if args.ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                                contig_name=contig,
                                                                region_start=pos - 1,
                                                                region_end=pos)
        if not pass_bed_region:
            truth_out_of_bed += 1
            del truth_variant_dict[key]
            continue

        if high_confident_only and key in low_qual_truth:
            continue

        pass_straed_region = len(strat_bed_tree_list) == 0 or sum([1 if is_region_in(tree=strat_bed_tree,
                                                                                     contig_name=contig,
                                                                                     region_start=pos - 1,
                                                                                     region_end=pos) else 0
                                                                   for strat_bed_tree in strat_bed_tree_list]) == len(
            strat_bed_tree_list)

        if not pass_straed_region and key in truth_variant_dict:
            del truth_variant_dict[key]
            truth_out_of_strat_bed += 1
            continue

        if key in low_af_truth:
            del truth_variant_dict[key]


    tp_snv, tp_ins, tp_del, fp_snv, fp_ins, fp_del, fn_snv, fn_ins, fn_del, fp_snv_truth, fp_ins_truth, fp_del_truth = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    truth_set = set()
    truth_snv, truth_ins, truth_del = 0, 0, 0
    query_snv, query_ins, query_del = 0, 0, 0
    pos_out_of_bed = 0

    fp_set = set()
    fn_set = set()
    fp_fn_set = set()
    tp_set = set()
    fp_qual_dict = defaultdict(float)
    tp_qual_dict = defaultdict(float)
    for key, vcf_infos in input_variant_dict.items():
        pos = key if args.ctg_name is not None else key[1]
        contig = args.ctg_name if args.ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                                contig_name=contig,
                                                                region_start=pos - 1,
                                                                region_end=pos)
        if not pass_bed_region:
            pos_out_of_bed += 1
            continue

        if high_confident_only and key in low_qual_truth:
            continue

        if key in low_af_truth:
            continue

        ref_base = vcf_infos.reference_bases
        alt_base = vcf_infos.alternate_bases[0]
        genotype = vcf_infos.genotype
        qual = vcf_infos.qual
        try:
            qual = float(qual) if qual is not None else qual
        except:
            qual = None
        is_snv = len(ref_base) == 1 and len(alt_base) == 1
        is_ins = len(ref_base) < len(alt_base)
        is_del = len(ref_base) > len(alt_base)

        if key not in truth_variant_dict and genotype != (0, 0):
            fp_snv = fp_snv + 1 if is_snv else fp_snv
            fp_ins = fp_ins + 1 if is_ins else fp_ins
            fp_del = fp_del + 1 if is_del else fp_del
            if is_snv:
                fp_set.add(key)
                fp_qual_dict[key] = qual
            if benchmark_indel and (is_ins or is_del):
                fp_set.add(key)
                fp_qual_dict[key] = qual

        if key in truth_variant_dict:
            vcf_infos = truth_variant_dict[key]
            truth_ref_base = vcf_infos.reference_bases
            truth_alt_base = vcf_infos.alternate_bases[0]
            truth_genotype = vcf_infos.genotype
            is_snv_truth = len(truth_ref_base) == 1 and len(truth_alt_base) == 1
            is_ins_truth = len(truth_ref_base) < len(truth_alt_base)
            is_del_truth = len(truth_ref_base) > len(truth_alt_base)

            if genotype == (0, 0) and truth_genotype == (0, 0):
                continue

            genotype_match = skip_genotyping or (truth_genotype == genotype)
            if truth_ref_base == ref_base and truth_alt_base == alt_base and genotype_match:
                tp_snv = tp_snv + 1 if is_snv else tp_snv
                tp_ins = tp_ins + 1 if is_ins else tp_ins
                tp_del = tp_del + 1 if is_del else tp_del
                if tp_snv or is_snv_truth:
                    tp_set.add(key)
                    tp_qual_dict[key] = qual
                if benchmark_indel and (is_ins or is_del):
                    tp_set.add(key)
                    tp_qual_dict[key] = qual

            else:
                fp_snv = fp_snv + 1 if is_snv else fp_snv
                fp_ins = fp_ins + 1 if is_ins else fp_ins
                fp_del = fp_del + 1 if is_del else fp_del

                fn_snv = fn_snv + 1 if is_snv_truth else fn_snv
                fn_ins = fn_ins + 1 if is_ins_truth else fn_ins
                fn_del = fn_del + 1 if is_del_truth else fn_del

                fp_set.add(key)
                fp_qual_dict[key] = qual
                fn_set.add(key)
                if is_snv or is_snv_truth:
                    fp_fn_set.add(key)
                if benchmark_indel and (is_ins_truth or is_del_truth):
                    fp_fn_set.add(key)

            truth_set.add(key)

    for key, vcf_infos in truth_variant_dict.items():
        pos = key if args.ctg_name is not None else key[1]
        contig = args.ctg_name if args.ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                                contig_name=contig,
                                                                region_start=pos - 1,
                                                                region_end=pos)

        if key in truth_set:
            continue
        if not pass_bed_region:
            continue

        if high_confident_only and key in low_qual_truth:
            continue

        if key in low_af_truth:
            continue

        truth_ref_base = vcf_infos.reference_bases
        truth_alt_base = vcf_infos.alternate_bases[0]
        truth_genotype = vcf_infos.genotype
        if truth_genotype == (0, 0):
            continue
        is_snv_truth = len(truth_ref_base) == 1 and len(truth_alt_base) == 1
        is_ins_truth = len(truth_ref_base) < len(truth_alt_base)
        is_del_truth = len(truth_ref_base) > len(truth_alt_base)

        fn_snv = fn_snv + 1 if is_snv_truth else fn_snv
        fn_ins = fn_ins + 1 if is_ins_truth else fn_ins
        fn_del = fn_del + 1 if is_del_truth else fn_del

        if is_snv_truth:
            fn_set.add(key)

        if benchmark_indel and (is_ins_truth or is_del_truth):
            fn_set.add(key)

    truth_indel = truth_ins + truth_del
    query_indel = query_ins + query_del
    tp_indel = tp_ins + tp_del
    fp_indel = fp_ins + fp_del
    fn_indel = fn_ins + fn_del

    snv_pre, snv_rec, snv_f1 = cal_metrics(tp=tp_snv, fp=fp_snv, fn=fn_snv)

    print("\n")
    print("[INFO] Total input records: {}, truth records: {}, records out of BED:{}".format(len(input_variant_dict),
                                                                                            len(truth_variant_dict),
                                                                                            pos_out_of_bed))
    add_tp_fn = False
    tp_fn = 'TP+FN' if add_tp_fn else ""
    tp_fn_count = tp_snv + fn_snv if add_tp_fn else ""
    print(''.join([item.ljust(13) for item in ["Type", 'Precision', 'Recall', "F1-score", 'TP', 'FP', 'FN', tp_fn]]),
          file=output_file)
    print(''.join([str(item).ljust(13) for item in ["SNV", snv_pre, snv_rec, snv_f1, tp_snv, fp_snv, fn_snv, tp_fn_count]]),
          file=output_file)
    if args.benchmark_indel:
        indel_pre, indel_rec, indel_f1 = cal_metrics(tp=tp_indel, fp=fp_indel, fn=fn_indel)
        ins_pre, ins_rec, ins_f1 = cal_metrics(tp=tp_ins, fp=fp_ins, fn=fn_ins)
        del_pre, del_rec, del_f1 = cal_metrics(tp=tp_del, fp=fp_del, fn=fn_del)
        print(''.join(
            [str(item).ljust(13) for item in ["INDEL", indel_pre, indel_rec, indel_f1, tp_indel, fp_indel, fn_indel, tp_indel+fn_indel]]),
              file=output_file)
        print(''.join(
            [str(item).ljust(13) for item in ["INS", ins_pre, ins_rec, ins_f1, tp_ins, fp_ins, fn_ins, tp_ins + fn_ins]]),
              file=output_file)
        print(''.join(
            [str(item).ljust(13) for item in ["DEL", del_pre, del_rec, del_f1, tp_del, fp_del, fn_del, tp_del + fn_del]]),
              file=output_file)

    if args.output_best_f1_score:
        results = output_best_cut_off(fp_qual_dict, tp_qual_dict, len(fn_set), use_int_cut_off=args.use_int_cut_off,
                                      add_tp_fn=add_tp_fn)
        best_match = results[0].copy()
        best_match[0] = 'SNV(Best F1)'
        print(
            ''.join(
                [str(item).ljust(13) if idx >= 4 or idx == 0 else ('%.4f' % item).ljust(13) for idx, item in
                 enumerate(best_match)]),
            file=output_file)

        if args.debug:
            print("")
            for result in results:
                print(''.join(
                    [str(item).ljust(13) if idx >= 4 or idx == 0 else ('%.4f' % item).ljust(13) for idx, item in
                     enumerate(result)]),
                    file=output_file)

    if args.roc_fn:
        if args.caller is None:
            fp_dict = dict([(key, float(input_variant_dict[key].qual)) for key in fp_set])
            tp_dict = dict([(key, float(input_variant_dict[key].qual)) for key in tp_set])
        elif args.caller.lower() == 'strelka2':
            fp_dict = {}
            for key in fp_set:
                somaticEVC = input_variant_dict[key].row_str.split('\t')[7].split(';')[-1]
                qual = float(somaticEVC.split('=')[1])
                fp_dict[key] = qual
            tp_dict = {}
            for key in tp_set:
                somaticEVC = input_variant_dict[key].row_str.split('\t')[7].split(';')[-1]
                qual = float(somaticEVC.split('=')[1])
                tp_dict[key] = qual

        elif args.caller.lower() == 'mutect2':
            fp_dict = {}
            for key in fp_set:
                TLOD = input_variant_dict[key].row_str.split('\t')[7].split(';')[-1]
                qual = float(TLOD.split('=')[1])
                fp_dict[key] = qual
            tp_dict = {}
            for key in tp_set:
                TLOD = input_variant_dict[key].row_str.split('\t')[7].split(';')[-1]
                qual = float(TLOD.split('=')[1])
                tp_dict[key] = qual

        elif args.caller.lower() == 'somaticsniper':
            fp_dict = {}
            for key in fp_set:
                SSC = input_variant_dict[key].row_str.split('\t')[10].split(':')[-1]
                qual = float(SSC)
                fp_dict[key] = qual
            tp_dict = {}
            for key in tp_set:
                SSC = input_variant_dict[key].row_str.split('\t')[10].split(':')[-1]
                qual = float(SSC)
                tp_dict[key] = qual

        elif args.caller.lower() == 'varnet':
            fp_dict = {}
            for key in fp_set:
                SCORE = input_variant_dict[key].row_str.split('\t')[7].split(';')[1]
                qual = float(SCORE.split('=')[1])
                fp_dict[key] = qual
            tp_dict = {}
            for key in tp_set:
                SCORE = input_variant_dict[key].row_str.split('\t')[7].split(';')[1]
                qual = float(SCORE.split('=')[1])
                tp_dict[key] = qual

        else:
            fp_dict = dict([(key, float(input_variant_dict[key].qual)) for key in fp_set])
            tp_dict = dict([(key, float(input_variant_dict[key].qual)) for key in tp_set])

        qual_list = sorted([float(qual) for qual in fp_dict.values()] + [qual for qual in tp_dict.values()],
                           reverse=True)

        tp_count = len(tp_set)
        roc_fn = open(args.roc_fn, 'w')
        for qual_cut_off in set(qual_list):
            pass_fp_count = sum([1 if float(qual) >= qual_cut_off else 0 for key, qual in fp_dict.items()])
            pass_tp_count = sum([1 if float(qual) >= qual_cut_off else 0 for key, qual in tp_dict.items()])
            fn_count = tp_count - pass_tp_count + fn_snv
            tmp_pre, tmp_rec, tmp_f1 = cal_metrics(tp=pass_tp_count, fp=pass_fp_count, fn=fn_count)
            roc_fn.write('\t'.join([str(round(item, 4)) for item in [qual_cut_off, tmp_pre, tmp_rec, tmp_f1]]) + '\n')
        roc_fn.close()

    if args.log_som is not None and os.path.exists(args.log_som):
        log_som = open(args.log_som)
        for row in log_som.readlines():
            if 'SNVs' not in row:
                continue

            columns = row.rstrip().split(',')
            tp, fp, fn, unk, ambi = [float(item) for item in columns[4:9]]
            recall, recall_lower, recall_upper, recall2 = [float(item) for item in columns[9:13]]
            precision, precision_lower, precision_upper = [float(item) for item in columns[13:16]]
            if int(tp_snv) != int(tp):
                print("True positives not match")
            if int(fp_snv) != int(fp):
                print("False positives not match")
            if int(fn_snv) != int(fn):
                print("False negatives not match")

    if output_dir is not None:
        if not os.path.exists(output_dir):
            subprocess.run("mkdir -p {}".format(output_dir), shell=True)
        candidate_types = ['fp', 'fn', 'fp_fn', 'tp']
        variant_sets = [fp_set, fn_set, fp_fn_set, tp_set]
        for vcf_type, variant_set in zip(candidate_types, variant_sets):
            vcf_fn = os.path.join(output_dir, '{}.vcf'.format(vcf_type))
            vcf_writer = VcfWriter(vcf_fn=vcf_fn, ctg_name=args.ctg_name, write_header=False)
            for key in variant_set:
                if key in input_variant_dict:
                    vcf_infos = input_variant_dict[key]
                elif key in truth_variant_dict:
                    vcf_infos = truth_variant_dict[key]
                else:
                    continue

                vcf_writer.write_row(row_str=vcf_infos.row_str)
            vcf_writer.close()
    if output_fn:
        output_file.close()


def main():
    parser = ArgumentParser(description="Compare input VCF with truth VCF")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="High confident BED region for benchmarking")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input VCF filename")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Truth VCF filename")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the contig to be proceed")

    parser.add_argument('--ctg_start', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctg_end', type=int, default=None,
                        help="The 1-based ending position of the sequence to be processed,")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output directory")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help="Filter variants with tag from the input VCF")

    parser.add_argument('--truth_filter_tag', type=str_none, default=None,
                        help="Filter variants with tag from the truth VCF")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument('--threads', type=int, default=8,
                        help="Max #threads to be used")

    # options for advanced users
    parser.add_argument('--min_af', type=float, default=None,
                        help="EXPERIMENTAL: Minimum VAF for a variant to be included in bechmarking")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum alt base count for a variant to be included in bechmarking")

    parser.add_argument('--min_coverage', type=int, default=4,
                        help="EXPERIMENTAL: Minimum coverage for a variant to be included in bechmarking")

    parser.add_argument('--strat_bed_fn', type=str, default=None,
                        help="EXPERIMENTAL: Genome stratifications v2 bed region")

    ## Output VCF filename
    parser.add_argument('--output_fn', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--low_af_path', type=str, default=None,
                        help=SUPPRESS)

    ## Only benchmark 'HighConf' tag in seqc VCF
    parser.add_argument('--high_confident_only', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--discard_fn_out_of_fp_bed', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--skip_genotyping', type=str2bool, default=True,
                        help="Skip calculating VCF genotype")

    parser.add_argument('--roc_fn', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--log_som', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--caller', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--output_best_f1_score', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--use_int_cut_off', type=str2bool, default=True,
                        help=SUPPRESS)

    parser.add_argument('--benchmark_indel', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--output_path', type=str, default=None,
                        help=SUPPRESS)

    ## output phase info
    parser.add_argument('--phase_output', type=str, default=None,
                        help=SUPPRESS)

    ## 0-> all, 1: phasable, 2 non-phasebale
    parser.add_argument('--validate_phase_only', type=str_none, default=None,
                        help=SUPPRESS)

    parser.add_argument('--min_qual', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--max_qual', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--naf_filter', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--debug', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    compare_vcf(args)


if __name__ == "__main__":
    main()
