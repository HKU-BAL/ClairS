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
import shlex
import logging
import random
import math

from subprocess import PIPE
from copy import deepcopy
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict, OrderedDict

import shared.param as param
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool, vcf_candidates_from
from shared.interval_tree import bed_tree_from, is_region_in

logging.basicConfig(format='%(message)s', level=logging.INFO)
BASES = set(list(BASE2NUM.keys()) + ["-"])
no_of_positions = param.no_of_positions
flanking_base_num = param.flankingBaseNum
channel_size = param.channel_size
BASE2NUMBER = dict(zip("ACGTURYSWKMBDHVN-", (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)))
NORMALIZE_NUM = param.NORMALIZE_NUM
ILMN_MAX_BQ = 40.0
ONT_MAX_BQ = 60.0
MAX_MQ = 60.0
MAX_AF = 1.0
STRAND_0 = -100
STRAND_1 = 100
HAP_TYPE = dict(zip((1, 0, 2), (30, 60, 90)))  # hap1 UNKNOWN H2 # should be better using
NORMAL_HAP_TYPE = dict(zip((1, 0, 2), (13, 25, 37)))  # set normal hap tag
# TUMOR_HAP_TYPE = dict(zip((1, 0, 2), (75, 88, 100)))  # set tumor hap tag
TUMOR_HAP_TYPE = dict(zip((1, 0, 2), (50, 75, 100)))  # set tumor hap tag

NORMAL_HAP_TYPE = dict(zip((1, 0, 2), (-90, -60, -30)))  # set normal hap tag
# TUMOR_HAP_TYPE = dict(zip((1, 0, 2), (75, 88, 100)))  # set tumor hap tag
TUMOR_HAP_TYPE = dict(zip((1, 0, 2), (30, 60, 90)))  # set tumor hap tag

# NORMAL_HAP_TYPE = dict(zip((1, 0, 2), (13, 50, 37)))  # set normal hap tag
# TUMOR_HAP_TYPE = dict(zip((1, 0, 2), (75, 100, 88)))  # set tumor hap tag
ACGT_NUM = dict(zip("ACGT+-*#N", (100, 25, 75, 50, -50, -100, -100, -100, 100)))


def add_gaussian_noise(data):
    noise_scale = 5

    def generate_normal_noise():
        u1 = random.uniform(0, 1)
        u2 = random.uniform(0, 1)
        z0 = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
        return z0 * noise_scale

    noise = round(generate_normal_noise())

    noisy_data = [max(1, min((item - noise), 90)) for item in data]

    return noisy_data


def normalize_bq(x, platform='ont'):
    MAX_BQ = ONT_MAX_BQ
    if platform == "ilmn":
        # only work for short reads
        x = x // 10 * 10
        MAX_BQ = ILMN_MAX_BQ
    return int(NORMALIZE_NUM * min(x, MAX_BQ) / MAX_BQ)


def normalize_mq(x):
    x = 0 if x <= 20 else (20 if x <= 40 else 60)
    return int(NORMALIZE_NUM * min(x, MAX_MQ) / MAX_MQ)


def get_chunk_id(candidates_bed_regions):
    chunk_id, chunk_num, bin_size = "", "", ""
    try:
        # path format: ctg_name:${chunk_id}_${index}_${bin_size}
        suffix = candidates_bed_regions.split('.')[-1]
        chunk_id, index, bin_size = suffix.split('_')
        index = str(int(index) + 1)
    except:
        return ""
    return "chunk {}-{}/{}".format(index, bin_size, chunk_id)


class Position(object):
    def __init__(self, pos, ref_base=None, alt_base=None, read_name_list=None, base_list=None, raw_base_quality=None,
                 raw_mapping_quality=None, af=None, depth=None, genotype=None, phase_set=None):
        self.pos = pos
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.read_name_list = read_name_list
        self.base_list = base_list
        self.raw_base_quality = raw_base_quality
        self.raw_mapping_quality = raw_mapping_quality
        self.af = af
        self.depth = depth
        self.read_channel = None
        self.mapping_quality = None
        self.update_info = False
        self.read_info = defaultdict()
        self.ref_seq = None
        self.alt_seq = None
        self.phase_set = phase_set
        self.genotype = genotype
        self.read_name_seq = defaultdict(str)

    def update_infos(self, is_tumor=False, mask_low_bq=False, platform='ont'):
        # only proceed when variant exists in candidate windows which greatly improves efficiency
        self.read_name_dict = dict(zip(self.read_name_list, self.base_list))
        self.update_info = True
        self.mapping_quality = [normalize_mq(phredscore2raw_score(item)) for item in self.raw_mapping_quality]
        # self.base_quality = [normalize_bq(phredscore2raw_score(item), platform) for item in self.raw_base_quality]
        self.base_quality = [normalize_bq(item, platform) for item in self.raw_base_quality]

        for read_name, base_info, bq, mq in zip(self.read_name_list, self.base_list, self.base_quality,
                                                self.mapping_quality):
            read_channel, ins_base, query_base = get_tensor_info(base_info, bq, self.ref_base, mask_low_bq, mq,
                                                                 is_tumor)
            self.read_info[read_name] = (read_channel, ins_base)


def phredscore2raw_score(qual):
    return ord(qual) - 33


def evc_base_from(base):
    if base == 'N':
        return 'A'
    elif base == 'n':
        return 'a'
    elif base in 'ACGTacgt':
        return base
    elif base.isupper():
        return 'A'
    else:
        return 'a'


def sorted_by_hap_read_name(center_pos, haplotag_dict, pileup_dict, hap_dict, max_depth, use_tensor_sample_mode=False):
    """
    Sort by reads haplotype after haplotag reads otherwise sort by read start position.
    center_pos: define the center candidate position for processing.
    haplotag_dict: dictionary (read name : hap type) which keep the read name and haplotype mapping.
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    hap_dict: similar to haplotag_dict, dictionary (pos: pos info) which keep the read name and haplotype mapping,
    while haplotype information directly acquire from BAM HP tag.
    platform: select maximum depth for each platform.
    """
    all_nearby_read_name = []
    start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num + 1
    for p in range(start_pos, end_pos):
        if p in pileup_dict.keys():
            all_nearby_read_name += pileup_dict[p].read_name_list
    all_nearby_read_name = list(OrderedDict.fromkeys(all_nearby_read_name))  # have sorted by order
    matrix_depth = max_depth
    if len(all_nearby_read_name) > matrix_depth and not use_tensor_sample_mode:
        # set same seed for reproducibility
        random.seed(0)
        indices = random.sample(range(len(all_nearby_read_name)), matrix_depth)
        all_nearby_read_name = [all_nearby_read_name[i] for i in sorted(indices)]
    sorted_read_name_list = []
    for order, read_name in enumerate(all_nearby_read_name):
        hap = max(haplotag_dict[read_name], hap_dict[read_name])  # no phasing is 0
        sorted_read_name_list.append((hap, order, read_name))

    sorted_read_name_list = sorted(sorted_read_name_list, key=lambda x: x[1])
    return sorted_read_name_list


def get_tensor_info(base_info, bq, ref_base, mask_low_bq=False, read_mq=None, is_tumor=False):
    """
    Create tensor information for each read level position.
    base_info: base information include all alternative bases.
    bq: normalized base quality.
    ref_base: reference_base: upper reference base for cigar calculation.
    read_mq: read mapping quality.
    """

    hap_type = TUMOR_HAP_TYPE[0] if is_tumor else NORMAL_HAP_TYPE[0]

    base, indel = base_info
    ins_base = ""
    query_base = ""
    read_channel = [0] * channel_size
    if base[0] in '*#':
        read_channel[1] = ACGT_NUM[base[0]]
        return read_channel, ins_base, query_base
    strand = STRAND_1
    if base[0] in 'ACGT':
        strand = STRAND_0
    ALT_BASE = 0

    base_upper = base.upper()
    if indel != '':
        ALT_BASE = ACGT_NUM[indel[0]]
    elif (base_upper != ref_base and base_upper in 'ACGT'):
        base_upper = evc_base_from(base_upper)
        ALT_BASE = ACGT_NUM[base_upper]
        if mask_low_bq and bq < 33 and ALT_BASE:  # bq < 20
            ALT_BASE = 0
            bq = 0

    REF_BASE = ACGT_NUM[ref_base]
    if len(indel) and indel[0] in '+-':
        if indel[0] == "+":
            ins_base = indel[1:].upper()
    read_channel[:6] = REF_BASE, ALT_BASE, strand, bq, read_mq, hap_type
    query_base = "" if base_upper not in "ACGT" else base_upper
    return read_channel, ins_base, query_base


def decode_pileup_bases(pos, pileup_bases, reference_base, minimum_snv_af_for_candidate, minimum_indel_af_for_candidate,
                        has_pileup_candidates, candidates_type_dict, is_tumor, platform="ont"):
    """
    Decode mpileup input string.
    pileup_bases: pileup base string for each position, include all mapping information.
    reference_base: upper reference base for cigar calculation.
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position.
    ref_seq: chunked reference sequence in window, start: center pos - flankingBaseNum, end: center + flankingBaseNum + 1.
    reference_sequence: reference sequence index by contig:start-end. 0-based.
    minimum_af_for_candidate: default minimum alleic frequency for candidate filtering, filter if below specific thredshold.
    has_pileup_candidates: if the candidate is directly obtained from pileup output, then no need to check the af filtering.
    """

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
    if has_pileup_candidates:
        if pos not in candidates_type_dict or not is_tumor:
            return base_list, None, True, 1.0
    pileup_dict = defaultdict(int)
    base_counter = Counter([''.join(item) for item in base_list])
    depth = 0
    for key, count in base_counter.items():
        if key[0].upper() in 'ACGT':
            pileup_dict[key[0].upper()] += count
            depth += count
        elif key[0] in "#*":
            depth += count
        if len(key) > 1 and key[1] == '+':
            pileup_dict['I'] += count
        elif len(key) > 1 and key[1] == '-':
            pileup_dict['D'] += count

    minimum_snv_af_for_candidate = minimum_snv_af_for_candidate if minimum_snv_af_for_candidate > 0 else param.min_af
    minimum_indel_af_for_candidate = minimum_indel_af_for_candidate if minimum_indel_af_for_candidate > 0 else param.min_af

    denominator = depth if depth > 0 else 1
    pileup_list = sorted(list(pileup_dict.items()), key=lambda x: x[1], reverse=True)

    pass_snv_af = False
    pass_indel_af = False

    for item, count in pileup_list:
        if pass_snv_af or pass_indel_af:
            break
        if item == reference_base:
            continue
        elif item[0] in 'ID':
            pass_indel_af = (pass_indel_af or (float(count) / denominator >= minimum_indel_af_for_candidate))
            continue
        pass_snv_af = pass_snv_af or (float(count) / denominator >= minimum_snv_af_for_candidate)

    af = (float(pileup_list[1][1]) / denominator) if len(pileup_list) > 1 else 0.0
    af = (float(pileup_list[0][1]) / denominator) if len(pileup_list) >= 1 and pileup_list[0][
        0] != reference_base else af

    pass_af = pass_snv_af or pass_indel_af

    return base_list, depth, pass_af, af


def find_tumor_alt_match(center_pos,
                         sorted_read_name_list,
                         pileup_dict,
                         truths_variant_dict,
                         proportion=None,
                         normal_output_bam_prefix='n',
                         tumor_output_bam_prefix='t',
                         use_real_sample=False):
    ref_base, alt_base = truths_variant_dict[center_pos].reference_bases, \
                         truths_variant_dict[center_pos].alternate_bases[0]
    is_ins = len(alt_base) > 1 and len(ref_base) == 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    is_snv = len(ref_base) == 1 and len(alt_base) == 1
    if use_real_sample:
        tumor_reads = [read_name for (hap, _, read_name) in sorted_read_name_list]
        matched_read_name_set = set()
        normal_read_name_set = None
        tumor_read_name_set = set(tumor_reads)
    else:
        if proportion is not None and float(proportion) == 1.0:
            # all reads are from tumor reads
            tumor_reads = [read_name for (hap, _, read_name) in sorted_read_name_list]
            normal_reads = []
        else:
            tumor_reads = [read_name for (hap, _, read_name) in sorted_read_name_list if read_name.startswith(tumor_output_bam_prefix)]
            normal_reads = [read_name for (hap, _, read_name) in sorted_read_name_list if read_name.startswith(normal_output_bam_prefix)]
        matched_read_name_set = set()
        normal_read_name_set = set(normal_reads)
        tumor_read_name_set = set(tumor_reads)
    for read_name in tumor_reads:
        if read_name in pileup_dict[center_pos].read_name_dict:
            base, indel = pileup_dict[center_pos].read_name_dict[read_name]
            base_upper = base.upper()
            if is_ins and base_upper + indel[1:].upper() == alt_base:
                matched_read_name_set.add(read_name)
            elif is_del and len(indel[1:].upper()) == len(ref_base[1:]):
                matched_read_name_set.add(read_name)
            elif is_snv and base_upper == alt_base:
                matched_read_name_set.add(read_name)
    return matched_read_name_set, normal_read_name_set, tumor_read_name_set


def generate_tensor(args,
                    ctg_name,
                    center_pos,
                    sorted_read_name_list,
                    pileup_dict,
                    ref_seq,
                    reference_sequence,
                    reference_start,
                    platform,
                    confident_bed_tree,
                    add_hetero_phasing,
                    is_tumor,
                    candidates_type_dict,
                    use_tensor_sample_mode=False,
                    truths_variant_dict=None,
                    proportion=None,
                    keep_phase_only=False,
                    hap_dict=None,
                    use_real_sample=False):
    """
    Generate full alignment input tensor
    ctg_name: provided contig name.
    center_pos: center position for full alignment generation, default window size = no_of_positions =
    flankingBaseNum + 1 + flankingBaseNum
    sorted_read_name_list: read name list which have been sorted by read start position and haplotype.
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    ref_seq: chunked reference sequence in window, start: center pos - flankingBaseNum, end: center + flankingBaseNum + 1.
    reference_sequence: reference sequence index by contig:start-end. 0-based.
    reference_base: upper reference base for cigar calculation.
    reference_start: upper reference base for cigar calculation.
    platform: platform for tensor shape, ont give a larger maximum depth compared with pb and illumina.
    confident_bed_tree: dictionary (contig name : intervaltree) for fast region query.
    add_no_phasing_data_training: boolean option to decide whether add no phasing data in training, we will
    resort the read and remove haplotype info when using this option.
    """

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    reference_base = ref_seq[flanking_base_num]
    tensor_depth = len(sorted_read_name_list)

    use_alt_base = False if platform != 'ilmn' else param.use_alt_base # need to update in indel

    if tensor_depth == 0:
        return None, None
    tensor = [[[0] * tensor_shape[2] for _ in range(tensor_shape[1])] for _ in range(tensor_depth)]
    start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num + 1
    insert_tuple = []

    # match deletion cases and bed format
    pass_confident_bed = not len(confident_bed_tree) or is_region_in(confident_bed_tree, ctg_name,
                                                                     center_pos - 2,
                                                                     center_pos + flanking_base_num + 1)
    if not pass_confident_bed:
        return None, None

    for p in range(start_pos, end_pos):
        if p in pileup_dict and not pileup_dict[p].update_info:
            pileup_dict[p].update_infos(is_tumor=is_tumor, mask_low_bq=args.mask_low_bq, platform=platform)
        for read_idx, read_name_info in enumerate(sorted_read_name_list):
            hap, read_order, read_name = read_name_info
            offset = p - start_pos
            if p in pileup_dict and read_name in pileup_dict[p].read_info:
                read_channel, ins_base = pileup_dict[p].read_info[read_name]
                tensor[read_idx][offset] = read_channel
                if ins_base != '' and p < end_pos - 1:
                    insert_tuple.append((read_idx, offset, ins_base, p))

    for read_idx, p, ins_base, center_p in insert_tuple:

        for ins_idx in range(min(len(ins_base), no_of_positions - p)):
            tensor[read_idx][ins_idx + p][6] = ACGT_NUM[ins_base[ins_idx]]

    matrix_depth = param.tumor_matrix_depth_dict[platform] if is_tumor else param.normal_matrix_depth_dict[platform]

    min_af_for_samping = 0.05
    max_af_for_sampling = 1.0
    chunk_read_size = 4
    min_tumor_support_read_num = param.min_tumor_support_read_num
    tensor_string_list = []
    alt_info_list = []
    # gradient = each reads
    if use_tensor_sample_mode:

        if len(sorted_read_name_list) > matrix_depth:
            # set same seed for reproducibility
            random.seed(0)
            sorted_read_name_list = random.sample(sorted_read_name_list, matrix_depth)

        tumor_reads_meet_alt_info_set, normal_read_name_set, tumor_read_name_set = find_tumor_alt_match(center_pos,
                                                                                                        sorted_read_name_list,
                                                                                                        pileup_dict,
                                                                                                        truths_variant_dict,
                                                                                                        proportion=proportion,
                                                                                                        normal_output_bam_prefix=args.normal_output_bam_prefix,
                                                                                                        tumor_output_bam_prefix=args.tumor_output_bam_prefix,
                                                                                                        use_real_sample=use_real_sample)
        if len(tumor_reads_meet_alt_info_set) == 0:
            print("No reads support tumor alternative in pos:{}".format(center_pos))
            return None, None
        tumor_read_name_list = list(tumor_read_name_set)
        normal_read_name_list = list(normal_read_name_set) if not use_real_sample else []
        tumor_reads_num = len(tumor_read_name_set)
        tumor_reads_meet_alt_info_num = len(tumor_reads_meet_alt_info_set)
        tumor_read_porportion = tumor_reads_meet_alt_info_num / float(tumor_reads_num)
        paired_reads_num = len(normal_read_name_set) if not use_real_sample else None
        sampled_reads_num_list = []

        if param.use_beta_subsampling:
            beta_acc_per = param.beta_acc_per
            sampled_reads_num_list.append(len(tumor_read_name_list))
            for s_idx in range(4):
                random.seed(center_pos + s_idx)
                random_af = random.random()
                af = None
                for pro_idx in range(len(beta_acc_per) - 1):
                    if random_af >= beta_acc_per[pro_idx] and random_af < beta_acc_per[pro_idx + 1]:
                        af = pro_idx / 100.0
                        break
                if af == None:
                    continue
                sampled_read_num = int((len(tumor_read_name_set) * af) / tumor_read_porportion)
                if sampled_read_num >= min_tumor_support_read_num and sampled_read_num <= len(tumor_read_name_set):
                    sampled_reads_num_list.append(sampled_read_num)
        else:
            for read_num in range(len(tumor_read_name_set)):
                if use_real_sample:
                    sampled_reads_num_list.append(read_num)
                else:
                    if read_num == 0 or paired_reads_num == 0 or read_num > matrix_depth:
                        continue
                    tumor_af = read_num / (read_num + paired_reads_num) * tumor_read_porportion
                    if tumor_af >= min_af_for_samping and tumor_af <= max_af_for_sampling:
                        sampled_reads_num_list.append(read_num)

            sampled_reads_num_list = sampled_reads_num_list[::chunk_read_size]

        sampled_reads_num_list = sorted(list(set(sampled_reads_num_list)))

        for read_num in sampled_reads_num_list:

            sampled_tumor_read_name_list = random.sample(tumor_read_name_list, read_num)
            sampled_tumor_read_name_meet_alt_set = set(sampled_tumor_read_name_list).intersection(
                tumor_reads_meet_alt_info_set)
            # in training, the candidates should have enough read support
            if len(sampled_tumor_read_name_meet_alt_set) < min_tumor_support_read_num:
                continue

            # check the strand info of sampled_tumor_read_name_list
            forward_strand = sum([1 for rn in sampled_tumor_read_name_meet_alt_set if rn.endswith('1')])
            if forward_strand == 0 or forward_strand == len(sampled_tumor_read_name_meet_alt_set):
                continue

            if use_alt_base:
                re_sorted_read_name_set = set(normal_read_name_list + tumor_read_name_list)
            else:
                re_sorted_read_name_set = set(normal_read_name_list + sampled_tumor_read_name_list)

            tmp_tensor = []
            tmp_depth = len(re_sorted_read_name_set)
            af_set = set()
            tmp_read_name_list = []
            alt_base_num = 0
            for row_idx, (hap, _, read_name) in enumerate(sorted_read_name_list):
                if read_name in re_sorted_read_name_set:
                    if use_alt_base and read_name not in sampled_tumor_read_name_meet_alt_set and truths_variant_dict[center_pos].alternate_bases[0] in ACGT_NUM.keys() and \
                            tensor[row_idx][flanking_base_num][1] == ACGT_NUM[truths_variant_dict[center_pos].alternate_bases[0]]:
                        tensor_copy = deepcopy(tensor[row_idx])
                        tensor_copy[flanking_base_num][1] = 0
                        tmp_tensor.append(tensor_copy)
                        alt_base_num += 1
                    else:
                        tmp_tensor.append(tensor[row_idx])
                    tmp_read_name_list.append(read_name)
            alt_dict = defaultdict(int)
            for read_name, (base, indel) in zip(pileup_dict[center_pos].read_name_list,
                                                pileup_dict[center_pos].base_list):
                if base in "#*":
                    continue
                if read_name not in re_sorted_read_name_set:
                    continue
                base_upper = base.upper()
                if indel != '':
                    if indel[0] == '+':
                        alt_dict['+' + base_upper + indel[1:].upper()] += 1
                    else:  # del
                        alt_dict[indel.upper()] += 1
                elif base.upper() != reference_base:
                    alt_dict[base.upper()] += 1

            if use_alt_base:
                alt_dict[truths_variant_dict[center_pos].alternate_bases[0]] -= alt_base_num
            alt_info = []
            af_infos = ' '.join([str(item) for item in sorted(list(af_set), reverse=True) if item != 0])
            for alt_type, alt_count in alt_dict.items():
                if alt_type[0] == '+':
                    alt_info.append(['I' + alt_type[1:].upper(), str(alt_count)])
                elif alt_type[0] == '-':
                    del_bases_num = len(alt_type[1:])
                    del_ref_bases = reference_sequence[
                                    center_pos - reference_start + 1:center_pos - reference_start + del_bases_num + 1]
                    alt_info.append(['D' + del_ref_bases, str(alt_count)])
                else:
                    alt_info.append(['X' + alt_type, str(alt_count)])
            alt_info = str(tmp_depth) + '-' + ' '.join(
                [' '.join([item[0], str(item[1])]) for item in alt_info]) + '-' + af_infos
            alt_info_list.append(alt_info)
            tensor_string_list.append(" ".join(
                (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in tmp_tensor))

            if add_hetero_phasing and ((candidates_type_dict[center_pos] != 'homo_somatic' or float(proportion) == 1.0) if not use_real_sample else True):
                HAP_TYPE = TUMOR_HAP_TYPE if is_tumor else NORMAL_HAP_TYPE
                unphased_num = TUMOR_HAP_TYPE[0] if is_tumor else NORMAL_HAP_TYPE[0]
                all_hap = [item[0] for item in sorted_read_name_list]
                # skip if no phased reads exist
                if sum(all_hap) != 0:

                    # require phasable haplotype for hetero somatic
                    if (candidates_type_dict[center_pos] == 'hetero_somatic' if not use_real_sample else True):
                        hap_counter = Counter([hap_dict[rn] for rn in sampled_tumor_read_name_meet_alt_set])
                        if hap_counter[1] * hap_counter[2] > 0 or (hap_counter[1] < 3 and hap_counter[2] < 3):
                            continue

                    sorted_phased_read_name_list = sorted(sorted_read_name_list, key=lambda x: (x[0], x[1]))
                    phase_read_name_index_mapping = [item[1] for item in sorted_phased_read_name_list]

                    phased_tensor = [deepcopy(tensor[read_idx]) for read_idx in phase_read_name_index_mapping]

                    for row_idx in range(len(phased_tensor)):
                        hap = sorted_phased_read_name_list[row_idx][0]
                        if hap in HAP_TYPE:
                            for p in range(no_of_positions):
                                if phased_tensor[row_idx][p][5] == unphased_num:
                                    phased_tensor[row_idx][p][5] = HAP_TYPE[hap]

                    phasing_tensor_string = " ".join(
                        (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in
                        phased_tensor)
                    if keep_phase_only:

                        tensor_string_list[-1] = phasing_tensor_string
                    else:
                        tensor_string_list.append(phasing_tensor_string)
                        alt_info_list.append(alt_info)

        variant_type = candidates_type_dict[center_pos] if center_pos in candidates_type_dict else 'unknown'

        return '\n'.join(["%s\t%d\t%s\t%s\t%s\t%s\t%s" % (
            ctg_name,
            center_pos,
            ref_seq,
            tensor_string,
            alt_info,
            "tumor" if is_tumor else "normal",
            variant_type
        ) for tensor_string, alt_info in zip(tensor_string_list, alt_info_list)]), ""

    else:
        alt_dict = defaultdict(int)
        depth, max_del_length = 0, 0
        for base, indel in pileup_dict[center_pos].base_list:
            if base in "#*":
                depth += 1
                continue
            depth += 1
            base_upper = base.upper()
            if indel != '':
                if indel[0] == '+':
                    alt_dict['+' + base_upper + indel[1:].upper()] += 1
                else:  # del
                    alt_dict[indel.upper()] += 1
                    max_del_length = max(len(indel), max_del_length)
            elif base.upper() != reference_base:
                alt_dict[base.upper()] += 1

        af_set = set()
        alt_info = []
        af_infos = ' '.join([str(item) for item in sorted(list(af_set), reverse=True) if item != 0])
        for alt_type, alt_count in alt_dict.items():
            if alt_type[0] == '+':
                alt_info.append(['I' + alt_type[1:].upper(), str(alt_count)])
            elif alt_type[0] == '-':
                del_bases_num = len(alt_type[1:])
                del_ref_bases = reference_sequence[
                                center_pos - reference_start + 1:center_pos - reference_start + del_bases_num + 1]
                alt_info.append(['D' + del_ref_bases, str(alt_count)])
            else:
                alt_info.append(['X' + alt_type, str(alt_count)])
        alt_info = str(depth) + '-' + ' '.join(
            [' '.join([item[0], str(item[1])]) for item in alt_info]) + '-' + af_infos
        tensor_string_list = [" ".join(
            (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in tensor)]

        if add_hetero_phasing and (candidates_type_dict[center_pos] != 'homo_somatic' if not use_real_sample else True):
            HAP_TYPE = TUMOR_HAP_TYPE if is_tumor else NORMAL_HAP_TYPE
            unphased_num = TUMOR_HAP_TYPE[0] if is_tumor else NORMAL_HAP_TYPE[0]
            all_hap = [item[0] for item in sorted_read_name_list]
            # skip if no phased reads exist
            if sum(all_hap) != 0:
                sorted_phased_read_name_list = sorted(sorted_read_name_list, key=lambda x: (x[0], x[1]))
                phase_read_name_index_mapping = [item[1] for item in sorted_phased_read_name_list]

                phased_tensor = [tensor[read_idx] for read_idx in phase_read_name_index_mapping]

                for row_idx in range(len(phased_tensor)):
                    hap = sorted_phased_read_name_list[row_idx][0]
                    if hap in HAP_TYPE:
                        for p in range(no_of_positions):
                            if phased_tensor[row_idx][p][5] == unphased_num:
                                phased_tensor[row_idx][p][5] = HAP_TYPE[hap]

                phasing_tensor_string = " ".join(
                    (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in
                    phased_tensor)

                if keep_phase_only:
                    tensor_string_list[-1] = phasing_tensor_string
                else:
                    tensor_string_list.append(phasing_tensor_string)

        variant_type = candidates_type_dict[center_pos] if center_pos in candidates_type_dict else 'unknown'

        return '\n'.join(["%s\t%d\t%s\t%s\t%s\t%s\t%s" % (
            ctg_name,
            center_pos,
            ref_seq,
            tensor_string,
            alt_info,
            "tumor" if is_tumor else "normal",
            variant_type
        ) for tensor_string in tensor_string_list]), alt_info


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def create_tensor(args):
    ctg_start = args.ctg_start
    ctg_end = args.ctg_end
    candidates_bed_regions = args.candidates_bed_regions
    fasta_file_path = args.ref_fn
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    bam_file_path = args.bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    tensor_can_output_path = args.tensor_can_fn
    is_candidates_bed_regions_given = candidates_bed_regions is not None
    add_hetero_phasing = phasing_info_in_bam = args.add_hetero_phasing
    keep_phase_only = args.keep_phase_only
    extend_bp = param.extend_bp
    add_bq_disturbance = args.add_bq_disturbance
    use_real_sample = args.use_real_sample
    minimum_snv_af_for_candidate = args.snv_min_af
    minimum_indel_af_for_candidate = args.indel_min_af
    min_coverage = args.min_coverage
    platform = args.platform
    confident_bed_fn = args.bed_fn
    alt_fn = args.alt_fn
    extend_bed = args.extend_bed
    is_extend_bed_file_given = extend_bed is not None
    min_mapping_quality = args.min_mq
    min_base_quality = args.min_bq
    vcf_fn = args.vcf_fn
    is_known_vcf_file_provided = vcf_fn is not None
    tensor_sample_mode = args.tensor_sample_mode
    is_tumor = args.is_tumor_sample
    global test_pos
    test_pos = None
    candidates_pos_set = set()
    candidates_type_dict = defaultdict(str)
    add_read_regions = True
    training_mode = args.training_mode
    truth_vcf_fn = args.truth_vcf_fn
    is_truth_vcf_provided = truth_vcf_fn is not None
    truths_variant_dict = {}
    proportion = args.proportion
    mask_low_bq = args.mask_low_bq

    if is_truth_vcf_provided:
        from shared.vcf import VcfReader
        if use_real_sample:
            unified_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, is_var_format=False,
                                           filter_tag="PASS")
        else:
            unified_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        truths_variant_dict = unified_vcf_reader.variant_dict

    if candidates_bed_regions:

        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (candidates_bed_regions)))
        candidate_file_path_output = candidate_file_path_process.stdout

        ctg_start, ctg_end = float('inf'), 0
        for row in candidate_file_path_output:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name: continue
            position = int(row[1]) + 1
            end = int(row[2]) + 1
            ctg_start = min(position, ctg_start)
            ctg_end = max(end, ctg_end)
            center = position + (end - position) // 2 - 1
            candidates_pos_set.add(center)
            variant_type = 'unknown'
            if len(row) == 4:
                variant_type = row[3]
            candidates_type_dict[center] = variant_type
        candidate_file_path_output.close()
        candidate_file_path_process.wait()

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    if is_known_vcf_file_provided:
        known_variants_list = vcf_candidates_from(vcf_fn=vcf_fn, contig_name=ctg_name)
        known_variants_set = set(known_variants_list)
    if not candidates_bed_regions and chunk_id is not None:

        """
        Whole genome calling option, acquire contig start end position from reference fasta index(.fai), then split the
        reference accroding to chunk id and total chunk numbers.
        """
        contig_length = 0
        with open(fai_fn, 'r') as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")

                contig_name = columns[0]
                if contig_name != ctg_name:
                    continue
                contig_length = int(columns[1])
        chunk_size = contig_length // chunk_num + 1 if contig_length % chunk_num else contig_length // chunk_num
        ctg_start = chunk_size * chunk_id  # 0-base to 1-base
        ctg_end = ctg_start + chunk_size

        # for illumina platform, the reads alignment is acquired after reads realignment from ReadsRealign.py
        if bam_file_path != "PIPE":
            bam_file_path += '.{}_{}'.format(ctg_start, ctg_end)
            add_read_regions = False
        if bam_file_path == "PIPE":
            add_read_regions = False

    # preparation for candidates near variants
    candidates_pos_set = set([item for item in candidates_pos_set if item >= ctg_start and item <= ctg_end])
    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    extend_start, extend_end = None, None
    if is_ctg_range_given:
        extend_start = ctg_start - no_of_positions
        extend_end = ctg_end + no_of_positions
        reads_regions.append(region_from(ctg_name=ctg_name, ctg_start=extend_start, ctg_end=extend_end))
        reference_start, reference_end = ctg_start - param.expandReferenceRegion, ctg_end + param.expandReferenceRegion
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=ctg_name, ctg_start=reference_start, ctg_end=reference_end))
    elif is_ctg_name_given:
        reads_regions.append(region_from(ctg_name=ctg_name))
        ref_regions.append(region_from(ctg_name=ctg_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools_execute_command,
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )
    if reference_sequence is None or len(reference_sequence) == 0:
        sys.exit("[ERROR] Failed to load reference sequence from file ({}).".format(fasta_file_path))

    phasing_option = " --output-extra HP" if phasing_info_in_bam else " "
    mq_option = ' --min-MQ {}'.format(min_mapping_quality)
    output_mq, output_read_name = True, True
    output_mq_option = ' --output-MQ ' if output_mq else ""
    output_read_name_option = ' --output-QNAME ' if output_read_name else ""
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    # pileup bed first
    bed_option = ' -l {}'.format(
        extend_bed) if is_extend_bed_file_given else ""
    bed_option = ' -l {}'.format(candidates_bed_regions) if is_candidates_bed_regions_given else bed_option
    flags_option = ' --excl-flags {} '.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    max_depth_option = ' --max-depth {}'.format(args.max_depth) if args.max_depth is not None else " "
    reads_regions_option = ' -r {}'.format(" ".join(reads_regions)) if add_read_regions else ""
    # print (add_read_regions, ctg_start, ctg_end, reference_start)
    stdin = None if bam_file_path != "PIPE" else sys.stdin
    bam_file_path = bam_file_path if bam_file_path != "PIPE" else "-"
    samtools_command = "{} mpileup {} --reverse-del".format(samtools_execute_command, bam_file_path) + \
                       output_read_name_option + output_mq_option + reads_regions_option + phasing_option + mq_option + bq_option + bed_option + flags_option + max_depth_option
    samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command), stdin=stdin)

    if tensor_can_output_path != "PIPE":
        tensor_can_fpo = open(tensor_can_output_path, "wb")
        tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(args.zstd)), stdin=PIPE, stdout=tensor_can_fpo)
    else:
        tensor_can_fp = TensorStdout(sys.stdout)

    if alt_fn:
        output_alt_fn = alt_fn
        alt_fpo = open(output_alt_fn, "wb")
        alt_fp = subprocess_popen(shlex.split("{} -c".format(args.zstd)), stdin=PIPE, stdout=alt_fpo)

    hap_dict = defaultdict(int)
    haplotag_dict = defaultdict(int)
    pileup_dict = defaultdict(str)
    extend_bp_distance = no_of_positions + param.extend_bp
    confident_bed_tree = bed_tree_from(bed_file_path=confident_bed_fn,
                                       contig_name=ctg_name,
                                       bed_ctg_start=extend_start,
                                       bed_ctg_end=extend_end)

    extend_bed_tree = bed_tree_from(bed_file_path=extend_bed,
                                    contig_name=ctg_name,
                                    bed_ctg_start=extend_start,
                                    bed_ctg_end=extend_end)

    def samtools_pileup_generator_from(samtools_mpileup_process):
        candidate_pos_list = sorted(list(candidates_pos_set))
        current_pos_index = 0
        has_pileup_candidates = len(candidates_pos_set)
        for row in samtools_mpileup_process.stdout:  # chr position N depth seq BQ read_name mapping_quality phasing_info
            columns = row.strip().split('\t')
            pos = int(columns[1])
            # pos that near bed region should include some indel cover in bed
            pass_extend_bed = not is_extend_bed_file_given or is_region_in(extend_bed_tree,
                                                                           ctg_name, pos - 1,
                                                                           pos + 1)
            pass_ctg_range = not ctg_start or (pos >= ctg_start and pos <= ctg_end)
            if not has_pileup_candidates and not pass_extend_bed and pass_ctg_range:
                continue
            pileup_bases = columns[4]
            raw_base_quality = columns[5]

            base_quality = [phredscore2raw_score(item) for item in raw_base_quality]

            if add_bq_disturbance:
                random.seed(ctg_start)
                new_base_quality = add_gaussian_noise(base_quality)
                base_quality = new_base_quality

            raw_mapping_quality = columns[6]
            read_name_list = columns[7].split(',')
            reference_base = reference_sequence[pos - reference_start].upper()
            if reference_base not in 'ACGT':
                continue
            base_list, depth, pass_af, af = decode_pileup_bases(pos=pos,
                                                                pileup_bases=pileup_bases,
                                                                reference_base=reference_base,
                                                                minimum_snv_af_for_candidate=minimum_snv_af_for_candidate,
                                                                minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
                                                                has_pileup_candidates=has_pileup_candidates,
                                                                candidates_type_dict=candidates_type_dict,
                                                                is_tumor=is_tumor)

            for b_idx, base in enumerate(base_list):
                if base[0] == '#' or (base[0] >= 'a' and base[0] <= 'z'):
                    read_name_list[b_idx] += '_1'  # reverse
                else:
                    read_name_list[b_idx] += '_0'  # forward

            if phasing_info_in_bam:
                phasing_info = columns[8].split(',')
                if len(read_name_list) != len(phasing_info):
                    continue
                else:
                    for hap_idx, hap in enumerate(phasing_info):
                        if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                            hap_dict[read_name_list[hap_idx]] = int(hap)

            if len(read_name_list) != len(base_list):
                continue

            if not is_known_vcf_file_provided and not has_pileup_candidates and reference_base in 'ACGT' and (
                    pass_af and depth >= min_coverage):
                candidate_pos_list.append(pos)

            if is_known_vcf_file_provided and not has_pileup_candidates and pos in known_variants_set:
                candidate_pos_list.append(pos)

            pileup_dict[pos] = Position(pos=pos,
                                        ref_base=reference_base,
                                        read_name_list=read_name_list,
                                        base_list=base_list,
                                        raw_base_quality=base_quality,
                                        raw_mapping_quality=raw_mapping_quality,
                                        af=af,
                                        depth=depth)

            # overlap_hetero_region = hetero_snv_tree.at(pos)

            if current_pos_index < len(candidate_pos_list) and pos - candidate_pos_list[
                current_pos_index] > extend_bp_distance:
                yield candidate_pos_list[current_pos_index]
                for pre_pos in sorted(pileup_dict.keys()):
                    if candidate_pos_list[current_pos_index] - pre_pos > extend_bp_distance:
                        del pileup_dict[pre_pos]
                    else:
                        break
                current_pos_index += 1
        while current_pos_index != len(candidate_pos_list):
            yield candidate_pos_list[current_pos_index]
            for pre_pos in sorted(pileup_dict.keys()):
                if candidate_pos_list[current_pos_index] - pre_pos > extend_bp_distance:
                    del pileup_dict[pre_pos]
                else:
                    break
            current_pos_index += 1

    samtools_pileup_generator = samtools_pileup_generator_from(samtools_mpileup_process)

    for pos in samtools_pileup_generator:
        if pos not in pileup_dict:
            continue

        use_tensor_sample_mode = is_tumor and tensor_sample_mode and (
                candidates_type_dict[pos] == 'homo_somatic' or candidates_type_dict[
            pos] == 'hetero_somatic') and pos in truths_variant_dict
        max_depth = param.tumor_matrix_depth_dict[platform] if is_tumor else param.normal_matrix_depth_dict[platform]
        sorted_read_name_list = sorted_by_hap_read_name(pos, haplotag_dict, pileup_dict, hap_dict, max_depth,
                                                        use_tensor_sample_mode)
        ref_seq = reference_sequence[
                  pos - reference_start - flanking_base_num: pos - reference_start + flanking_base_num + 1].upper()

        tensor, alt_info = generate_tensor(args=args,
                                           ctg_name=ctg_name,
                                           center_pos=pos,
                                           sorted_read_name_list=sorted_read_name_list,
                                           pileup_dict=pileup_dict,
                                           ref_seq=ref_seq,
                                           reference_sequence=reference_sequence,
                                           reference_start=reference_start,
                                           platform=platform,
                                           confident_bed_tree=confident_bed_tree,
                                           add_hetero_phasing=add_hetero_phasing,
                                           is_tumor=is_tumor,
                                           candidates_type_dict=candidates_type_dict,
                                           use_tensor_sample_mode=use_tensor_sample_mode,
                                           truths_variant_dict=truths_variant_dict,
                                           proportion=proportion if not use_real_sample else None,
                                           keep_phase_only=keep_phase_only,
                                           hap_dict=hap_dict,
                                           use_real_sample=use_real_sample)
        if not tensor:
            continue

        tensor_can_fp.stdin.write(tensor)
        tensor_can_fp.stdin.write("\n")
        if alt_fn:
            alt_fp.stdin.write('\t'.join([ctg_name + ' ' + str(pos), alt_info]) + '\n')

    samtools_mpileup_process.stdout.close()
    samtools_mpileup_process.wait()

    if tensor_can_output_path != "PIPE":
        tensor_can_fp.stdin.close()
        tensor_can_fp.wait()
        tensor_can_fpo.close()

    if alt_fn:
        alt_fp.stdin.close()
        alt_fp.wait()
        alt_fpo.close()


def main():
    parser = ArgumentParser(description="Generate variant candidate full-alignment tensors for training")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input, default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--snv_min_af', type=float, default=param.snv_min_af,
                        help="Minimum snv allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctg_start', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--ctg_end', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctg_name and/or (--ctg_start, --ctg_end) are set")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for advanced users
    parser.add_argument('--min_coverage', type=float, default=param.min_coverage,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=None,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--alt_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")


    # options for internal process control
    ## Minimum indel allele frequency for a site to be considered as a candidate site
    parser.add_argument('--indel_min_af', type=float, default=1.0,
                        help=SUPPRESS)

    ## Path to the 'zstd' compression
    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help=SUPPRESS)

    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=str2bool, default=0,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--candidates_bed_regions', type=str, default=None,
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform. Greatly boost indel performance in trade of running time
    parser.add_argument('--add_hetero_phasing', type=str2bool, default=0,
                        help=SUPPRESS)

    ## only keep phasing tensor for model training
    parser.add_argument('--keep_phase_only', type=str2bool, default=0,
                        help=SUPPRESS)

    ## apply tensor sample mode in training
    parser.add_argument('--tensor_sample_mode', type=str2bool, default=1,
                        help=SUPPRESS)

    ## apply tensor sample mode in training
    parser.add_argument('--is_tumor_sample', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--normal_output_bam_prefix', type=str, default='n',
                        help="Normal output BAM prefix")

    parser.add_argument('--tumor_output_bam_prefix', type=str, default='t',
                        help="Tumor output BAM prefix")

    parser.add_argument('--mask_low_bq', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--training_mode', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--add_phasing_info', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--proportion', type=float, default=1.0,
                        help=SUPPRESS)

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--add_bq_disturbance', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--use_real_sample', type=str2bool, default=0,
                        help=SUPPRESS)

    args = parser.parse_args()

    create_tensor(args)


if __name__ == "__main__":
    main()
