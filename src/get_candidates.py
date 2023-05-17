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

import shlex
import os
import sys
import math
import argparse
import random
import shlex
import math
import subprocess

from os.path import isfile, abspath
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
from subprocess import PIPE
from shared.vcf import VcfWriter
from shared.bed import BedWriter
from collections import defaultdict
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict

from shared.interval_tree import bed_tree_from, is_region_in
import shared.param as param
from shared.utils import AltInfos, str2bool


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def vcf_reader(vcf_fn, contig_name, bed_tree=None, add_hetero_pos=False):
    homo_variant_set = set()
    variant_set = set()
    variant_info = defaultdict()
    homo_variant_info = defaultdict()
    hetero_variant_set = set()
    hetero_variant_info = defaultdict()

    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    for row in unzip_process.stdout:
        row = row.rstrip()
        if row[0] == '#':
            continue
        columns = row.strip().split('\t')

        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        variant_set.add(pos)
        variant_info[pos] = (ref_base, alt_base)
        if bed_tree and not is_region_in(tree=bed_tree, contig_name=contig_name, region_start=pos):
            continue

        genotype_info = columns[9].split(':')
        genotype = genotype_info[0]
        g1, g2 = genotype.replace('|', "/").split('/')
        if g1 == "1" and g2 == "1":
            homo_variant_set.add(pos)
            homo_variant_info[pos] = (ref_base, alt_base)
        if add_hetero_pos and (g1 == "1" and g2 == "0" or g1 == "0" and g2 == "1"):
            hetero_variant_set.add(pos)
            hetero_variant_info[pos] = (ref_base, alt_base)

    return homo_variant_set, homo_variant_info, hetero_variant_set, hetero_variant_info, variant_set, variant_info


def get_ref_candidates(fn, contig_name=None, bed_tree=None, variant_info=None):
    ref_cans_dict = defaultdict(AltInfos)
    if os.path.exists(fn):
        fn_list = [fn]
    else:
        fn = fn.split('/')
        directry, file_prefix = '/'.join(fn[:-1]), fn[-1]
        fn_list = [os.path.join(directry, f) for f in os.listdir(directry) if f.startswith(file_prefix)]

    if len(fn_list) == 0:
        print('[ERROR] No file prefix')
        return ref_cans_dict
    for fn in fn_list:
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (fn)))
        for row in unzip_process.stdout:
            if row[0] == '#':
                continue
            columns = row.rstrip().split('\t')
            ctg_name, pos = columns[0], int(columns[1])
            if contig_name and ctg_name != contig_name:
                continue
            if not is_region_in(bed_tree, contig_name, pos):
                continue

            ref_base, depth, af_infos, alt_infos = columns[2:6]
            tumor_infos = columns[6] if len(columns) > 6 else ""

            af_list = af_infos.split(',')
            alt_dict = dict([[item.split(':')[0], float(item.split(':')[1])] for item in alt_infos.split(' ')])
            tumor_alt_dict = dict(
                [[item.split(':')[0], float(item.split(':')[1])] for item in tumor_infos.split(' ')]) if len(
                tumor_infos) else dict()
            ref_cans_dict[pos] = AltInfos(pos=pos,
                                          ref_base=ref_base,
                                          depth=depth,
                                          af_list=af_list,
                                          alt_dict=alt_dict,
                                          tumor_alt_dict=tumor_alt_dict)
        unzip_process.stdout.close()
        unzip_process.wait()
    return ref_cans_dict


def find_most_frequent_candidate(args, alt_info_dict, ref_base):
    alt_list = sorted(list([k, v] for k, v in alt_info_dict.items() if k != ref_base), key=lambda x: x[1], reverse=True)

    for alt_base, af in alt_list:
        if alt_base == ref_base:
            continue
        if not args.select_indel_candidates and len(alt_base) > 1:
            continue
        return alt_base, af

    return None, None


def find_candidate_match(alt_info_dict, ref_base, alt_base):
    # alt_base = alt_base[0]
    alt_list = sorted(list(alt_info_dict.items()), key=lambda x: x[1], reverse=True)
    max_af = alt_list[0][1] if len(alt_info_dict) > 0 else None
    # snp
    if len(ref_base) == 1 and len(alt_base) == 1:
        if alt_base in alt_info_dict:
            return alt_info_dict[alt_base], 'snp', max_af
    # insertion
    if len(ref_base) == 1 and len(alt_base) > 1:
        ab = alt_base[0] + '+' + alt_base[1:]
        if ab in alt_info_dict:
            return alt_info_dict[ab], 'ins', max_af
    # deletion
    if len(ref_base) > 1 and len(alt_base) == 1:
        ab = ref_base[0] + '-' + 'N' * len(ref_base[1:])
        if ab in alt_info_dict:
            return alt_info_dict[ab], 'del', max_af

    if len(alt_info_dict) > 0:
        return None, 'signal', max_af
    return None, "", None


def filter_germline_candidates(args, truths, variant_info, alt_dict, paired_alt_dict, gen_vcf, INFO, platform='ont'):
    # alt_dict: tumor
    # pair alt_dict: normal
    filtered_truths = []
    low_confident_germline_truths = []
    truth_not_pass_af = 0
    germline_filtered_by_af_distance = 0
    germline_af_gap = 0.1 if platform == 'ont' else 0.1
    min_germline_af = 0.1 if platform == 'ont' else None
    select_indel_candidates = args.select_indel_candidates
    for pos, variant_type in truths:
        if pos not in alt_dict:
            truth_not_pass_af += 1
            if gen_vcf:
                low_confident_germline_truths.append((pos, variant_type + INFO + 'germline_no_read_support'))

            continue

        # discard the germline variants, if normal and tumor have large distinct
        if pos in paired_alt_dict:
            ref_base, alt_base = variant_info[pos]
            if select_indel_candidates:
                if len(ref_base) == 1 and len(alt_base) == 1:
                    continue

            pair_af, vt, max_af = find_candidate_match(alt_info_dict=paired_alt_dict[pos].alt_dict, ref_base=ref_base,
                                                       alt_base=alt_base)
            paired_alt_dict[pos].max_candidate_af = max_af
            paired_alt_dict[pos].support_alternative_af = pair_af

            af, vt, max_af = find_candidate_match(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base,
                                                  alt_base=alt_base)
            alt_dict[pos].max_candidate_af = max_af
            alt_dict[pos].support_alternative_af = af
            # pair_af: af in normal, af: af in tumor
            if pair_af is None or af is None or math.fabs(pair_af - af) > germline_af_gap or (
                    min_germline_af is not None and (af < min_germline_af or pair_af < min_germline_af)):
                germline_filtered_by_af_distance += 1
                if gen_vcf:
                    low_confident_germline_truths.append((pos, variant_type + INFO + 'germline_af_gap'))
                continue

            if pos in paired_alt_dict and paired_alt_dict[pos].support_alternative_af is None:
                paired_alt_dict[pos].support_alternative_af = pair_af
            if pos in alt_dict and alt_dict[pos].support_alternative_af is None:
                alt_dict[pos].support_alternative_af = af
            filtered_truths.append([pos, variant_type])

    print(
        '[INFO] {} truth germline variants filtered by no read support in normal BAM: {}, filtered by large AF distance between normal and tumor: {}'.format(
            INFO, truth_not_pass_af, germline_filtered_by_af_distance))

    return filtered_truths, low_confident_germline_truths


def filter_reference_candidates(args, truths, alt_dict, paired_alt_dict, gen_vcf, INFO):
    # we filter the reference candidates that have small gap between normal and tumor pair
    filtered_truths = []
    low_confident_reference_candidates = []
    truth_not_pass_af = 0
    reference_filtered_by_af_distance = 0
    for pos, variant_type in truths:
        if pos not in alt_dict:
            truth_not_pass_af += 1
            if gen_vcf:
                low_confident_reference_candidates.append((pos, variant_type + INFO + 'reference_no_read_support'))
            continue

        tumor_af, pair_af = None, None
        ref_base = alt_dict[pos].ref_base
        most_frequent_alt_base, tumor_af = find_most_frequent_candidate(args=args,
                                                                        alt_info_dict=alt_dict[pos].alt_dict,
                                                                        ref_base=ref_base)

        if most_frequent_alt_base is None or tumor_af is None:
            continue
        tumor_af = float(tumor_af)
        if pos not in paired_alt_dict:
            filtered_truths.append([pos, variant_type])
        else:
            pair_af, vt, pair_max_af = find_candidate_match(alt_info_dict=paired_alt_dict[pos].alt_dict,
                                                            ref_base=ref_base, alt_base=most_frequent_alt_base)
            if pair_af is None:
                filtered_truths.append([pos, variant_type])
            elif pair_af >= 0.05 or (pair_af > tumor_af * 0.2):
                reference_filtered_by_af_distance += 1
                if gen_vcf:
                    low_confident_reference_candidates.append((pos, variant_type + INFO + 'reference_af_gap'))
                continue

            filtered_truths.append([pos, variant_type])

        if pos in paired_alt_dict and paired_alt_dict[pos].support_alternative_af is None:
            paired_alt_dict[pos].support_alternative_af = pair_af
        if pos in alt_dict and alt_dict[pos].support_alternative_af is None:
            alt_dict[pos].support_alternative_af = tumor_af

    print(
        '[INFO] {} reference calls filtered by no read support in normal BAM: {}, filtered by large AF distance between normal and tumor: {}'.format(
            INFO, truth_not_pass_af, reference_filtered_by_af_distance))

    return filtered_truths, low_confident_reference_candidates


def filter_somatic_candidates(args, truths, variant_info, alt_dict, paired_alt_dict, gen_vcf=False, INFO="Homo"):
    filtered_truths = []
    truth_not_pass_af = 0
    truth_filter_in_normal = 0
    truth_filter_with_low_af = 0
    truth_filter_with_low_coverage = 0
    min_af_for_tumor = 0.03
    max_af_in_normal = 0.03
    min_tumor_support_read_num = param.min_tumor_support_read_num
    # af_gap_for_errors = 0.15
    low_confident_truths = []
    skip_no_read_support = True
    skip_high_af_in_normal = True
    skip_low_af_in_tumor = True
    skip_low_coverage_in_tumor = True
    for pos, variant_type in truths:
        if pos not in alt_dict:
            truth_not_pass_af += 1
            if gen_vcf:
                low_confident_truths.append((pos, variant_type + INFO + '_no_read_support'))
            if skip_no_read_support:
                continue

        normal_af, tumor_reads_af = None, None

        if pos in paired_alt_dict:
            ref_base, alt_base = variant_info[pos]
            if args.select_indel_candidates:
                if len(ref_base) == 1 and len(alt_base) == 1:
                    continue

            # very high af with same alt_base in normal is not suitable as candidate for training
            normal_af, vt, max_af = find_candidate_match(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base,
                                                         alt_base=alt_base)
            if (normal_af is not None and normal_af > max_af_in_normal) or (
                    max_af is not None and max_af > max_af_in_normal):
                truth_filter_in_normal += 1
                if gen_vcf:
                    low_confident_truths.append((pos, variant_type + INFO + '_high_af_in_normal'))
                if skip_high_af_in_normal:
                    continue

        tumor_alt_dict = alt_dict[pos].tumor_alt_dict if alt_dict[pos].tumor_alt_dict is not None else alt_dict[
            pos].alt_dict
        if len(tumor_alt_dict):
            ref_base, alt_base = variant_info[pos]
            tumor_reads_af, vt, max_af = find_candidate_match(alt_info_dict=tumor_alt_dict, ref_base=ref_base,
                                                              alt_base=alt_base)
            if tumor_reads_af is None or tumor_reads_af < min_af_for_tumor:
                truth_filter_with_low_af += 1
                if gen_vcf:
                    low_confident_truths.append((pos, variant_type + INFO + '_low_af_in_tumor'))
                if skip_low_af_in_tumor:
                    continue
            tumor_coverage = int(alt_dict[pos].depth) * tumor_reads_af
            if tumor_coverage < min_tumor_support_read_num:
                truth_filter_with_low_coverage += 1
                if gen_vcf:
                    low_confident_truths.append((pos, variant_type + INFO + '_low_coverage_in_tumor'))
                if skip_low_coverage_in_tumor:
                    continue

        if pos in paired_alt_dict and paired_alt_dict[pos].support_alternative_af is None:
            paired_alt_dict[pos].support_alternative_af = normal_af
        if pos in alt_dict and alt_dict[pos].support_alternative_af is None:
            alt_dict[pos].support_alternative_af = tumor_reads_af

        filtered_truths.append([pos, variant_type])

    print(
        '[INFO] {} truth variants filtered by high AF in normal BAM: {}, filtered by low AF {}:{}, filtered by no/low read support:{}/{}'.format(
            INFO, truth_filter_in_normal, min_af_for_tumor, truth_filter_with_low_af, truth_not_pass_af,
            truth_filter_with_low_coverage))
    return filtered_truths, low_confident_truths


def get_candidates(args):
    contig_name = args.ctg_name
    bed_fn = args.bed_fn
    bed_tree = bed_tree_from(bed_fn, contig_name)
    normal_vcf_fn = args.normal_vcf_fn
    tumor_vcf_fn = args.tumor_vcf_fn
    maximum_non_variant_ratio = args.maximum_non_variant_ratio
    normal_reference_cans_fn = args.normal_reference_cans
    tumor_reference_cans_fn = args.tumor_reference_cans
    add_hetero_pos = args.add_hetero_pos
    split_bed_size = param.split_bed_size
    flanking_base_num = args.flanking_base_num if args.flanking_base_num else param.flankingBaseNum
    split_folder = args.split_folder
    output_vcf_fn = args.output_vcf_fn
    gen_vcf = output_vcf_fn is not None
    sample_normal_af = args.sample_normal_af
    ref_fn = args.ref_fn
    platform = args.platform
    proportion = args.proportion
    synthetic_coverage = args.synthetic_coverage
    output_bed_fn = args.output_bed_fn
    output_fp_bed_regions = output_bed_fn is not None
    exclude_flanking_truth = args.exclude_flanking_truth

    normal_homo_variant_set, normal_homo_variant_info, normal_hetero_variant_set, normal_hetero_variant_info, normal_variant_set, normal_variant_info = vcf_reader(
        vcf_fn=normal_vcf_fn, contig_name=contig_name, bed_tree=bed_tree, add_hetero_pos=add_hetero_pos)
    tumor_homo_variant_set, tumor_homo_variant_info, tumor_hetero_variant_set, tumor_hetero_variant_info, tumor_variant_set, tumor_variant_info = vcf_reader(
        vcf_fn=tumor_vcf_fn, contig_name=contig_name, bed_tree=bed_tree, add_hetero_pos=add_hetero_pos)
    tumor_alt_dict = get_ref_candidates(fn=tumor_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree,
                                        variant_info=tumor_variant_info)
    normal_alt_dict = get_ref_candidates(fn=normal_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree,
                                         variant_info=tumor_variant_info)

    normal_ref_cans_list = [pos for pos in normal_alt_dict if
                            pos not in tumor_variant_set and pos not in normal_variant_set]
    tumor_ref_cans_list = [pos for pos in tumor_alt_dict if
                           pos not in tumor_variant_set and pos not in normal_variant_set]
    intersection_pos_set = normal_homo_variant_set.intersection(tumor_homo_variant_set)

    same_alt_pos_set = set()
    for pos in intersection_pos_set:
        if normal_homo_variant_info[pos] != tumor_homo_variant_info[pos]:
            continue
        same_alt_pos_set.add(pos)

    hetero_list_with_same_repre = []
    if add_hetero_pos:
        for pos, (ref_base, alt_base) in tumor_hetero_variant_info.items():
            if pos in normal_hetero_variant_set:
                ref_base_2, alt_base_2 = normal_hetero_variant_info[pos]
                if ref_base == ref_base_2 and alt_base == alt_base_2:
                    hetero_list_with_same_repre.append(pos)

    homo_germline = [(item, 'homo_germline') for item in list(same_alt_pos_set)]
    hetero_germline = [(item, 'hetero_germline') for item in hetero_list_with_same_repre]
    ref_list = tumor_ref_cans_list
    if sample_normal_af is not None:
        random.seed(0)
        ref_list += random.sample(normal_ref_cans_list, int(len(normal_ref_cans_list) * sample_normal_af))

    references = [(item, 'ref') for item in ref_list]
    homo_germline, homo_low_confident_germline_truths = filter_germline_candidates(args=args,
                                                                                   truths=homo_germline,
                                                                                   variant_info=tumor_variant_info,
                                                                                   alt_dict=tumor_alt_dict,
                                                                                   paired_alt_dict=normal_alt_dict,
                                                                                   gen_vcf=gen_vcf,
                                                                                   INFO="Homo",
                                                                                   platform=platform)
    # need add hetero, otherwise, in real case, the performance of hetero variants are too bad
    hetero_germline, hetero_low_confident_germline_truths = filter_germline_candidates(args=args,
                                                                                       truths=hetero_germline,
                                                                                       variant_info=tumor_variant_info,
                                                                                       alt_dict=tumor_alt_dict,
                                                                                       paired_alt_dict=normal_alt_dict,
                                                                                       gen_vcf=gen_vcf,
                                                                                       INFO="Hetero",
                                                                                       platform=platform)

    add_germline = True
    if not add_germline:
        # exclude germline variants if we exclude germline variants into training
        # exclude gerlmine variants when create VCF and fp BED region
        homo_germline = []
        hetero_germline = []

    references, low_confident_reference_candidates = filter_reference_candidates(args=args,
                                                                                 truths=references,
                                                                                 alt_dict=tumor_alt_dict,
                                                                                 paired_alt_dict=normal_alt_dict,
                                                                                 gen_vcf=gen_vcf,
                                                                                 INFO="Ref")

    if maximum_non_variant_ratio is not None:
        # random sample reference calls in training mode, with seed
        random.seed(0)
        references = random.sample(references, int(len(references) * maximum_non_variant_ratio))
        random.seed(0)
        homo_germline = random.sample(homo_germline, int(len(homo_germline) * maximum_non_variant_ratio))
        random.seed(0)
        hetero_germline = random.sample(hetero_germline, int(len(hetero_germline) * maximum_non_variant_ratio))

    if args.use_reference_candidates_only:
        homo_germline = []
        hetero_germline = []

    fp_list = homo_germline + references + hetero_germline

    homo_somatic_set = sorted(list(tumor_homo_variant_set - normal_variant_set))
    hetero_somatic_set = sorted(list(tumor_hetero_variant_set - normal_variant_set))
    homo_somatic = [(item, 'homo_somatic') for item in homo_somatic_set]
    # skip hetero variant here
    hetero_somatic = [(item, 'hetero_somatic') for item in hetero_somatic_set] if add_hetero_pos else []

    if args.use_reference_candidates_only:
        homo_somatic = []
        hetero_somatic = []

    homo_somatic, homo_low_confident_truths = filter_somatic_candidates(args=args,
                                                                        truths=homo_somatic,
                                                                        variant_info=tumor_variant_info,
                                                                        alt_dict=tumor_alt_dict,
                                                                        paired_alt_dict=normal_alt_dict,
                                                                        gen_vcf=gen_vcf)

    hetero_somatic, hetero_low_confident_truths = filter_somatic_candidates(args=args,
                                                                            truths=hetero_somatic,
                                                                            variant_info=tumor_variant_info,
                                                                            alt_dict=tumor_alt_dict,
                                                                            paired_alt_dict=normal_alt_dict,
                                                                            gen_vcf=gen_vcf,
                                                                            INFO="Hetero")

    # exclude nearby truth in training
    if exclude_flanking_truth:
        all_truth_pos_list = sorted([item[0] for item in homo_somatic + hetero_somatic + fp_list])

        def exclude_flanking(pos_list):
            exclude_truth_set = set()
            truth_pos_list = sorted([item[0] for item in pos_list])
            for p_idx, pos in enumerate(truth_pos_list):
                start = pos - flanking_base_num
                end = pos + flanking_base_num + 1
                for idx in range(len(all_truth_pos_list)):
                    if all_truth_pos_list[idx] == pos:
                        continue
                    if all_truth_pos_list[idx] < start:
                        continue
                    elif all_truth_pos_list[idx] >= end:
                        break
                    random.seed(pos)
                    if random.random() > 0.7:
                        continue
                    exclude_truth_set.add(pos)
            return exclude_truth_set

        homo_exclude_truth_set = exclude_flanking(homo_somatic)
        hetero_exclude_truth_set = exclude_flanking(hetero_somatic)
        homo_somatic = [item for item in homo_somatic if item[0] not in homo_exclude_truth_set]
        hetero_somatic = [item for item in hetero_somatic if item[0] not in hetero_exclude_truth_set]

    tp_list = homo_somatic + hetero_somatic

    pos_list = sorted(fp_list + tp_list, key=lambda x: x[0])

    if gen_vcf:
        vcf_writer = VcfWriter(vcf_fn=output_vcf_fn, ref_fn=ref_fn, ctg_name=contig_name, show_ref_calls=True)
        for pos, variant_type in pos_list + homo_low_confident_truths + hetero_low_confident_truths + homo_low_confident_germline_truths + hetero_low_confident_germline_truths:
            genotype = '1/1' if (variant_type == 'homo_somatic' or variant_type == 'hetero_somatic') else '0/0'
            filter_tag = "PASS" if genotype == '1/1' else variant_type
            if genotype == "1/1":
                ref_base, alt_base = tumor_variant_info[pos]
            elif pos in tumor_alt_dict:
                ref_base = alt_base = tumor_alt_dict[pos].ref_base
            else:
                alt_base = ref_base = normal_alt_dict[pos].ref_base if pos in normal_alt_dict else (
                    tumor_alt_dict[pos] if pos in tumor_alt_dict else "N")

            normal_coverage = int(normal_alt_dict[pos].depth) if pos in normal_alt_dict else 0
            tumor_coverage = int(tumor_alt_dict[pos].depth) if pos in tumor_alt_dict else 0
            normal_af = normal_alt_dict[pos].support_alternative_af if pos in normal_alt_dict else 0
            tumor_af = tumor_alt_dict[pos].support_alternative_af if pos in tumor_alt_dict else 0

            normal_coverage = -1 if normal_coverage is None else normal_coverage
            tumor_coverage = -1 if tumor_coverage is None else tumor_coverage
            normal_af = -1 if normal_af is None else normal_af
            tumor_af = -1 if tumor_af is None else tumor_af
            vcf_writer.write_row(POS=pos,
                                 REF=ref_base,
                                 ALT=alt_base,
                                 QUAL=100,
                                 FILTER=filter_tag,
                                 GT=genotype,
                                 DP=tumor_coverage,
                                 AF=tumor_af,
                                 NAF=normal_af,
                                 NDP=normal_coverage,
                                 TDP=tumor_coverage,
                                 VT=variant_type)
        vcf_writer.close()

    if output_fp_bed_regions:
        bed_writer = BedWriter(bed_fn=output_bed_fn)
        for pos, _ in pos_list:
            bed_writer.write_row(ctg_name=contig_name,
                                 start_pos=pos - 1,
                                 end_pos=pos)
        bed_writer.close()

    all_candidates_regions = []
    region_num = len(pos_list) // split_bed_size + 1 if len(
        pos_list) % split_bed_size else len(pos_list) // split_bed_size

    for idx in range(region_num):
        # a windows region for create tensor # samtools mpileup not include last position
        split_output = pos_list[idx * split_bed_size: (idx + 1) * split_bed_size]

        split_output = [(item[0] - flanking_base_num, item[0] + flanking_base_num + 2, item[1]) for item in
                        split_output]

        output_path = os.path.join(split_folder,
                                   '{}.{}_{}_por{}_cov{}'.format(contig_name, idx, region_num, int(proportion * 100),
                                                                 int(synthetic_coverage)))
        all_candidates_regions.append(output_path + ' ' + str(proportion) + ' ' + str(synthetic_coverage))
        with open(output_path, 'w') as output_file:
            output_file.write('\n'.join(
                ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), x[2]]) for x in
                 split_output]) + '\n')  # bed format

    all_candidate_path = os.path.join(split_folder,
                                      'CANDIDATES_FILE_{}_{}_{}'.format(contig_name, proportion, synthetic_coverage))
    with open(all_candidate_path, 'w') as output_file:
        output_file.write('\n'.join(all_candidates_regions) + '\n')

    extra_info = "homo/hetetro somatic exclude by flanking position:{}/{}".format(len(homo_exclude_truth_set), len(
        hetero_exclude_truth_set)) if exclude_flanking_truth else ""
    print(
        '[INFO] {}-{} homo germline:{} hetero germline:{} references:{} homo somatic:{} hetero somatic:{} {}\n'.format(
            contig_name, proportion, len(homo_germline), len(hetero_germline), len(references), len(homo_somatic),
            len(hetero_somatic), extra_info))


def main():
    parser = ArgumentParser(description="Get pair cofident candidates for model training")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--normal_vcf_fn', type=str, default=None,
                        help="Normal candidate sites VCF file input")

    parser.add_argument('--tumor_vcf_fn', type=str, default=None,
                        help="Tumor candidate sites VCF file input")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--normal_reference_cans', type=str, default=None,
                        help="Normal reference cans information")

    parser.add_argument('--tumor_reference_cans', type=str, default=None,
                        help="Tumor reference cans information")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--output_dir', type=str, default="",
                        help="Output directory")

    parser.add_argument('--select_indel_candidates', type=str2bool, default=0,
                        help="Get Indel candidates instead of SNV candidates")

    # options for debug purpose
    parser.add_argument('--sample_normal_af', type=float, default=None,
                        help="DEBUG: Sample a proportion of normal candidates")

    parser.add_argument('--proportion', type=float, default=None,
                        help="Synthetic proportion of input BAM")

    parser.add_argument('--synthetic_coverage', type=int, default=None,
                        help="Synthetic coverage in training")

    # options for debug purpose
    parser.add_argument('--add_hetero_pos', type=str2bool, default=0,
                        help="Add hetero candidates into training")

    parser.add_argument('--exclude_flanking_truth', type=str2bool, default=1,
                        help="Exclude truths in a flanking window into training")

    parser.add_argument('--use_reference_candidates_only', type=str2bool, default=0,
                        help="Exclude truths in a flanking window into training")

    ## Output VCF path
    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help=SUPPRESS)
    ## Output BED path
    parser.add_argument('--output_bed_fn', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--split_folder', type=str, default=None,
                        help=SUPPRESS)

    # Flanking window size
    parser.add_argument('--flanking_base_num', type=int, default=None,
                        help=SUPPRESS)

    ## Maximum non-variant ratio against variant in the training data
    parser.add_argument('--maximum_non_variant_ratio', type=float, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    get_candidates(args)


if __name__ == "__main__":
    main()
