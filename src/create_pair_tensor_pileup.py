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
import heapq
from subprocess import PIPE
from itertools import product
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict

import shared.param as param
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool, vcf_candidates_from
from shared.interval_tree import bed_tree_from, is_region_in
from src.create_tensor import get_chunk_id

logging.basicConfig(format='%(message)s', level=logging.INFO)
BASES = set(list(BASE2NUM.keys()) + ["-"])
no_of_positions = param.no_of_positions
flanking_base_num = param.flankingBaseNum


BASE2NUMBER = dict(zip("ACGTURYSWKMBDHVN-", (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)))
NORMALIZE_NUM = param.NORMALIZE_NUM
HAP_TYPE = dict(zip((1, 0, 2), (30, 60, 90)))  # hap1 UNKNOWN H2
ACGT_NUM = dict(zip("ACGT+-*#N", (100, 25, 75, 50, -50, -100, 0, 0, 100)))

#           0    1    2    3    4    5    6    7     8    9    10   11  12   13    14  15   16    17  18      19      20      21
channel = ['A', 'C', 'G', 'T', 'I', 'I1', 'D', 'D1', '*', 'a', 'c', 'g','t', 'i', 'i1','d', 'd1','#']
channel += [
 'ALMQ', 'CLMQ', 'GLMQ', 'TLMQ', 'aLMQ', 'cLMQ', 'gLMQ', 'tLMQ', 'ALBQ', 'CLBQ', 'GLBQ', 'TLBQ', 'aLBQ', 'cLBQ', 'gLBQ', 'tLBQ']

phase_channel = [
 'AHP1', 'CHP1', 'GHP1', 'THP1', 'aHP1', 'cHP1', 'gHP1', 'tHP1', 'AHP2', 'CHP2', 'GHP2', 'THP2', 'aHP2', 'cHP2', 'gHP2', 'tHP2']
channel_size = len(channel)
BASE2INDEX = dict(zip(channel, tuple(range(channel_size))))
from src.create_tensor_pileup import base_index, phase_channel


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


def decode_pileup_bases(args,
                        pos,
                        pileup_bases,
                        reference_base,
                        minimum_snp_af_for_candidate,
                        minimum_indel_af_for_candidate,
                        has_pileup_candidates,
                        candidates_type_dict,
                        is_tumor,
                        mapping_quality,
                        base_quality,
                        phasing_info=None,
                        chunk_ref_seq=None,
                        platform="ont"):
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
    pileup_tensor = [0] * (channel_size if phasing_info is None else (channel_size + len(phase_channel)))
    is_candidate = pos in candidates_type_dict
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

    pileup_dict = defaultdict(int)
    base_counter = Counter([''.join(item) for item, mq in zip(base_list, mapping_quality) if mq >= 20])
    low_mq_base_counter = Counter([''.join(item) for item, mq in zip(base_list, mapping_quality) if mq < 20])
    low_bq_base_counter = Counter([''.join(item) for item, bq in zip(base_list, base_quality) if bq < (30 if platform == 'ont' else 10)])
    if phasing_info is not None:
        for b, hap in zip(base_list, phasing_info):
            base = ''.join(b)
            if hap in '12' and base in 'ACGTacgt':
                pileup_tensor[channel_size + base_index[base + "HP" + hap]] += 1

    depth, max_ins_0, max_del_0, max_ins_1, max_del_1 = 0, 0, 0, 0, 0
    max_del_length = 0
    alt_info_set = set()
    forward_alt_info_dict = defaultdict(int)
    reverse_alt_info_dict = defaultdict(int)
    forward_ref_count = 0
    reverse_ref_count = 0
    forward_depth = 0
    reverse_depth = 0
    for key, count in base_counter.items():
        if len(key) == 1:
            if key.upper() in 'ACGT':
                pileup_dict[key.upper()] += count
                if is_candidate and key.upper() != reference_base:
                    alt_info_set.add('X' + key.upper())
                    if key in "ACGT":
                        forward_alt_info_dict['X' + key.upper()] += count
                    else:
                        reverse_alt_info_dict['X' + key.upper()] += count
                if is_candidate and key.upper() == reference_base:
                    if key in 'ACGTN':
                        forward_ref_count += count
                    else:
                        reverse_ref_count = count
                pileup_tensor[BASE2INDEX[key]] += count
            elif key in '#*':
                pileup_tensor[BASE2INDEX[key]] += count

            if key in 'ACGTN*':
                forward_depth += count
            else:
                reverse_depth += count

        elif key[1] == '+':
            indel_len = len(key[2:])
            if indel_len > args.max_indel_length:
                continue
            if key[0] in "ACGTN":
                forward_depth += count
            else:
                reverse_depth += count

            pileup_dict['I'] += count
            if is_candidate:
                k = 'I' + key[0].upper() + key[2:].upper()
                alt_info_set.add(k)
                if key[0] in 'ACGTN':
                    forward_alt_info_dict[k] += count
                else:
                    reverse_alt_info_dict[k] += count
            # two strand
            if key[0] in 'ACGTN*':
                pileup_tensor[BASE2INDEX["I"]] += count
                max_ins_0 = max(max_ins_0, count)
            else:
                pileup_tensor[BASE2INDEX["i"]] += count
                max_ins_1 = max(max_ins_1, count)
        elif key[1] == '-':
            indel_len = len(key[1:])
            if indel_len > args.max_indel_length:
                continue

            if key[0] in "N*ACGT":
                forward_depth += count
            else:
                reverse_depth += count

            pileup_dict['D'] += count
            if is_candidate:
                info = chunk_ref_seq[:len(key[1:])]
                alt_info_set.add('D' + info)
                if key[0] in 'N*ACGT':
                    forward_alt_info_dict['D' + info] += count
                else:
                    reverse_alt_info_dict['D' + info] += count
            max_del_length = max(max_del_length, len(key[1:]))
            # two strand
            if key[0] in 'N*ACGT':
                pileup_tensor[BASE2INDEX["D"]] += count
                max_del_0 = max(max_del_0, count)
            else:
                pileup_tensor[BASE2INDEX["d"]] += count
                max_del_1 = max(max_del_1, count)
    if is_candidate and forward_ref_count + reverse_ref_count > 0:
        forward_alt_info_dict['R' + reference_base] = forward_ref_count
        reverse_alt_info_dict['R' + reference_base] = reverse_ref_count

    alt_info = str(forward_depth) + '/' + str(reverse_depth) + '-' + ' '.join(
                [' '.join([k, str(forward_alt_info_dict[k]) + '/' + str(reverse_alt_info_dict[k])]) for k in alt_info_set]) + '-'
    pileup_tensor[BASE2INDEX['I1']] = max_ins_0
    pileup_tensor[BASE2INDEX['i1']] = max_ins_1
    pileup_tensor[BASE2INDEX['D1']] = max_del_0
    pileup_tensor[BASE2INDEX['d1']] = max_del_1

    for key, count in low_mq_base_counter.items():
        if key.upper() in 'ACGT':
            pileup_tensor[BASE2INDEX[key+'LMQ']] += count

    for key, count in low_bq_base_counter.items():
        if key.upper() in 'ACGT':
            pileup_tensor[BASE2INDEX[key+'LBQ']] += count

    pileup_tensor[BASE2INDEX[reference_base]] = sum([pileup_tensor[BASE2INDEX[item]] for item in 'ACGT'])
    pileup_tensor[BASE2INDEX[reference_base.lower()]] = sum([pileup_tensor[BASE2INDEX[item]] for item in 'acgt'])
    pileup_tensor[BASE2INDEX[reference_base+'LMQ']] = sum([pileup_tensor[BASE2INDEX[item+'LMQ']] for item in 'ACGT'])
    pileup_tensor[BASE2INDEX[reference_base.lower()+'LMQ']] = sum([pileup_tensor[BASE2INDEX[item+'LMQ']] for item in 'acgt'])
    pileup_tensor[BASE2INDEX[reference_base+'LBQ']] = sum([pileup_tensor[BASE2INDEX[item+'LBQ']] for item in 'ACGT'])
    pileup_tensor[BASE2INDEX[reference_base.lower()+'LBQ']] = sum([pileup_tensor[BASE2INDEX[item+'LBQ']] for item in 'acgt'])


    if has_pileup_candidates:
    #     if pos not in candidates_type_dict or not is_tumor:
        return pileup_tensor, base_list, None, True, 1.0, alt_info

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

    minimum_snp_af_for_candidate = minimum_snp_af_for_candidate if minimum_snp_af_for_candidate > 0 else param.min_af
    minimum_indel_af_for_candidate = minimum_indel_af_for_candidate if minimum_indel_af_for_candidate > 0 else param.min_af_dict[platform]

    denominator = depth if depth > 0 else 1
    pileup_list = sorted(list(pileup_dict.items()), key=lambda x: x[1], reverse=True)

    pass_snp_af = False
    pass_indel_af = False

    for item, count in pileup_list:
        if pass_snp_af or pass_indel_af:
            break
        if item == reference_base:
            continue
        elif item[0] in 'ID':
            pass_indel_af = (pass_indel_af or (float(count) / denominator >= minimum_indel_af_for_candidate))
            continue
        pass_snp_af = pass_snp_af or (float(count) / denominator >= minimum_snp_af_for_candidate)

    af = (float(pileup_list[1][1]) / denominator) if len(pileup_list) > 1 else 0.0
    af = (float(pileup_list[0][1]) / denominator) if len(pileup_list) >= 1 and pileup_list[0][
        0] != reference_base else af

    pass_af = pass_snp_af or pass_indel_af

    return pileup_tensor, base_list, depth, pass_af, af, alt_info



class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def update_hetero_ref(pos, reference_sequence, reference_start, extend_bp, alt_base):
    # if need phasing option enables, will store reference sequence near hetero snp candidate.
    ref_start = pos - extend_bp
    ref_end = pos + extend_bp + 1
    ref_seq = reference_sequence[ref_start - reference_start: ref_end - reference_start]
    alt_seq = ref_seq[:extend_bp] + alt_base + ref_seq[extend_bp + 1:]
    return ref_seq, alt_seq

def get_normal_set(alt_fn):
    normal_pos_set = set()
    file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (alt_fn)))
    file_path_output = file_path_process.stdout
    for row in file_path_output:
        ctg_name, pos = row.rstrip().split(maxsplit=2)
        normal_pos_set.add(int(pos))
    file_path_output.close()
    file_path_process.wait()
    return normal_pos_set

def heapq_merge_generator_from(normal_bam_pileup_generator, tumor_bam_pileup_generator, skip_if_normal_empty=True):
    normal_candidates_set = set()
    tumor_candidates_set = set()
    for pos, is_tumor in heapq.merge(normal_bam_pileup_generator, tumor_bam_pileup_generator):
        if is_tumor:
            if pos in normal_candidates_set or not skip_if_normal_empty:
                yield pos
                normal_candidates_set.discard(pos)
            else:
                continue
        else:
            normal_candidates_set.add(pos)

    for pos in tumor_candidates_set:
        if pos in normal_candidates_set:
            yield pos

def get_key_list(input_dict, normal_tensor_infos_dict, tumor_tensor_infos_dict, shuffle = True):

    normal_index_list = range(len(input_dict['normal'][0]))
    tumor_index_list = range(len(input_dict['tumor'][0]))
    for normal_idx, tumor_idx in product(normal_index_list, tumor_index_list):
        normal_tensor, normal_alt_info = normal_tensor_infos_dict[0][normal_idx], normal_tensor_infos_dict[1][normal_idx]
        tumor_tensor, tumor_alt_info = tumor_tensor_infos_dict[0][normal_idx], tumor_tensor_infos_dict[1][normal_idx]
        yield (normal_tensor, normal_alt_info, tumor_tensor, tumor_alt_info)


def create_tensor(args):
    ctg_start = args.ctg_start
    ctg_end = args.ctg_end
    candidates_bed_regions = args.candidates_bed_regions
    fasta_file_path = args.ref_fn
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    normal_bam_file_path = args.normal_bam_fn
    tumor_bam_file_path = args.tumor_bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    tensor_can_output_path = args.tensor_can_fn
    is_candidates_bed_regions_given = candidates_bed_regions is not None
    minimum_snp_af_for_candidate = args.snv_min_af
    minimum_indel_af_for_candidate = args.indel_min_af
    min_coverage = args.min_coverage
    platform = args.platform
    confident_bed_fn = file_path_from(args.bed_fn, allow_none=True, exit_on_not_found=False)
    extend_bed = file_path_from(args.extend_bed, allow_none=True, exit_on_not_found=False)
    is_extend_bed_file_given = extend_bed is not None
    min_mapping_quality = args.min_mq
    min_base_quality = args.min_bq if args.min_bq is not None else param.min_bq_dict[platform]
    vcf_fn = args.vcf_fn
    is_known_vcf_file_provided = vcf_fn is not None
    phasing_info_in_bam = args.phase_tumor and args.platform == 'ont'
    args.max_indel_length = param.max_indel_length if args.max_indel_length is None else args.max_indel_length

    candidates_pos_set = set()
    candidates_type_dict = defaultdict(str)
    add_read_regions = True
    flanking_base_num = param.flankingBaseNum if args.flanking is None else args.flanking
    no_of_positions = 2 * flanking_base_num + 1

    truth_vcf_fn = args.truth_vcf_fn
    is_truth_vcf_provided = truth_vcf_fn is not None
    truths_variant_dict = {}
    if is_truth_vcf_provided:
        from shared.vcf import VcfReader
        unified_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        truths_variant_dict = unified_vcf_reader.variant_dict

    if candidates_bed_regions:

        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (candidates_bed_regions)))
        candidate_file_path_output = candidate_file_path_process.stdout

        ctg_start, ctg_end = float('inf'), 0
        for row in candidate_file_path_output:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name:
                continue
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

    normal_phasing_option = " "
    tumor_phasing_option = " --output-extra HP " if phasing_info_in_bam else " "
    samtools_view_min_mq = 0
    # mq_option = ' --min-MQ {}'.format(min_mapping_quality)
    mq_option = ' --min-MQ {}'.format(samtools_view_min_mq)
    output_mq, output_read_name = True, False
    output_mq_option = ' --output-MQ ' if output_mq else ""
    output_read_name_option = ' --output-QNAME ' if output_read_name else ""
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    # pileup bed first
    bed_option = ' -l {}'.format(
        extend_bed) if is_extend_bed_file_given else ""
    bed_option = ' -l {}'.format(candidates_bed_regions) if is_candidates_bed_regions_given else bed_option
    flags_option = ' --excl-flags {}'.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    max_depth_option = ' --max-depth {}'.format(args.max_depth) if args.max_depth is not None else ""

    reads_regions_option = ' -r {}'.format(" ".join(reads_regions)) if add_read_regions else ""
    # print (add_read_regions, ctg_start, ctg_end, reference_start)

    samtools_command = "{} mpileup --reverse-del".format(samtools_execute_command) + \
                       output_read_name_option + output_mq_option + reads_regions_option + mq_option + bq_option + bed_option + flags_option + max_depth_option
    samtools_mpileup_normal_process = subprocess_popen(
        shlex.split(samtools_command + normal_phasing_option + ' ' + normal_bam_file_path), stderr=PIPE)

    samtools_mpileup_tumor_process = subprocess_popen(
        shlex.split(samtools_command + tumor_phasing_option + ' ' + tumor_bam_file_path), stderr=PIPE)


    if tensor_can_output_path != "PIPE":
        tensor_can_fpo = open(tensor_can_output_path, "wb")
        tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(args.zstd)), stdin=PIPE, stdout=tensor_can_fpo)
    else:
        tensor_can_fp = TensorStdout(sys.stdout)

    extend_bp_distance = no_of_positions + param.extend_bp

    extend_bed_tree = bed_tree_from(bed_file_path=extend_bed,
                                    contig_name=ctg_name,
                                    bed_ctg_start=extend_start,
                                    bed_ctg_end=extend_end)

    normal_pileup_tensors = [[0] * channel_size] * (extend_end - extend_start + flanking_base_num)
    tumor_pileup_tensors = [[0] * channel_size] * (extend_end - extend_start + flanking_base_num)

    normal_alt_info_dict = defaultdict()
    tumor_alt_info_dict = defaultdict()

    def samtools_pileup_generator_from(samtools_mpileup_process, is_tumor=True):
        candidate_pos_list = sorted(list(candidates_pos_set))
        current_pos_index = 0
        has_pileup_candidates = len(candidates_pos_set)
        alt_info_dict = tumor_alt_info_dict if is_tumor else normal_alt_info_dict
        pileup_tensors = tumor_pileup_tensors if is_tumor else normal_pileup_tensors

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
            raw_mapping_quality = columns[6]
            reference_base = evc_base_from(reference_sequence[pos - reference_start]).upper()
            if reference_base not in 'ACGT':
                continue

            mapping_quality = [ord(mq) - 33 for mq in raw_mapping_quality]
            base_quality = [ord(mq) - 33 for mq in raw_base_quality]

            if phasing_info_in_bam and is_tumor:
                phasing_info = columns[7].split(',')
            else:
                phasing_info = None

            chunk_ref_seq = reference_sequence[pos - reference_start: pos - reference_start + args.max_indel_length].upper()

            pileup_tensor, base_list, depth, pass_af, af, alt_info = decode_pileup_bases(args=args,
                                                                                         pos=pos,
                                                                                         pileup_bases=pileup_bases,
                                                                                         reference_base=reference_base,
                                                                                         minimum_snp_af_for_candidate=minimum_snp_af_for_candidate,
                                                                                         minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
                                                                                         has_pileup_candidates=has_pileup_candidates,
                                                                                         candidates_type_dict=candidates_type_dict,
                                                                                         mapping_quality=mapping_quality,
                                                                                         base_quality=base_quality,
                                                                                         phasing_info=phasing_info,
                                                                                         chunk_ref_seq=chunk_ref_seq,
                                                                                         is_tumor=is_tumor)

            offset = pos - extend_start
            pileup_tensors[offset] = pileup_tensor
            if pos in candidates_type_dict:
                alt_info_dict[pos] = alt_info

            if not is_known_vcf_file_provided and not has_pileup_candidates and reference_base in 'ACGT' and (
                    pass_af and depth >= min_coverage):
                candidate_pos_list.append(pos)

            if is_known_vcf_file_provided and not has_pileup_candidates and pos in known_variants_set:
                candidate_pos_list.append(pos)

            if current_pos_index < len(candidate_pos_list) and pos - candidate_pos_list[
                current_pos_index] > extend_bp_distance:
                yield (candidate_pos_list[current_pos_index], is_tumor)

                current_pos_index += 1
        while current_pos_index != len(candidate_pos_list):
            yield (candidate_pos_list[current_pos_index], is_tumor)
            current_pos_index += 1

    normal_bam_pileup_generator = samtools_pileup_generator_from(samtools_mpileup_process=samtools_mpileup_normal_process,is_tumor=False)
    tumor_bam_pileup_generator = samtools_pileup_generator_from(samtools_mpileup_process=samtools_mpileup_tumor_process)

    tensor_count = 0
    for pos in heapq_merge_generator_from(normal_bam_pileup_generator=normal_bam_pileup_generator, tumor_bam_pileup_generator=tumor_bam_pileup_generator):
        ref_seq = reference_sequence[
                  pos - reference_start - flanking_base_num: pos - reference_start + flanking_base_num + 1].upper()
        start_index = pos - flanking_base_num - extend_start
        end_index = start_index + no_of_positions
        if start_index < 0 or end_index >= extend_end - extend_start:
            continue
        if sum([1 for tensor in normal_pileup_tensors[start_index:end_index] if len(tensor) == 0]) > 0:
            continue

        if sum([1 for tensor in tumor_pileup_tensors[start_index:end_index] if len(tensor) == 0]) > 0:
            continue

        variant_type = candidates_type_dict[pos] if pos in candidates_type_dict else 'unknown'

        tensor_infos_dict = defaultdict()
        normal_tensor_string_list = [" ".join(" ".join("%d" % x for x in innerlist) for innerlist in normal_pileup_tensors[start_index:end_index])]
        tumor_tensor_string_list = [" ".join(" ".join("%d" % x for x in innerlist) for innerlist in tumor_pileup_tensors[start_index:end_index])]
        if pos not in normal_alt_info_dict or pos not in tumor_alt_info_dict:
            continue
        tensor_infos_dict['normal'] = (normal_tensor_string_list, [normal_alt_info_dict[pos]])
        tensor_infos_dict['tumor'] = (tumor_tensor_string_list, [tumor_alt_info_dict[pos]])

        for tensor_infos in get_key_list(tensor_infos_dict, tensor_infos_dict['normal'], tensor_infos_dict['tumor']):
            normal_tensor_string, normal_alt_info, tumor_tensor_string, tumor_alt_info = tensor_infos

            tensor = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                ctg_name,
                pos,
                ref_seq,
                normal_tensor_string,
                normal_alt_info,
                tumor_tensor_string,
                tumor_alt_info,
                variant_type)
            tensor_can_fp.stdin.write(tensor)
            tensor_count += 1
    samtools_mpileup_normal_process.stdout.close()
    samtools_mpileup_normal_process.wait()
    samtools_mpileup_tumor_process.stdout.close()
    samtools_mpileup_tumor_process.wait()
    if tensor_can_output_path != "PIPE":
        tensor_can_fp.stdin.close()
        tensor_can_fp.wait()
        tensor_can_fpo.close()

    chunk_info = get_chunk_id(candidates_bed_regions)
    print("[INFO] {} {} Tensors generated: {}".format(ctg_name, chunk_info, tensor_count))


def main():
    parser = ArgumentParser(description="Generate tumor-normal pair pileup tensors for calling")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input, required")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--snv_min_af', type=float, default=param.snv_min_af,
                        help="Minimum snp allele frequency for a site to be considered as a candidate site, default: %(default)f")

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

    parser.add_argument('--min_bq', type=int, default=None,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=None,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--alt_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    # options for internal process control
    ## Path to the 'zstd' compression
    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help=SUPPRESS)

    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=str2bool, default=0,
                        help=SUPPRESS)

    ## Minimum indel allele frequency for a site to be considered as a candidate site
    parser.add_argument('--indel_min_af', type=float, default=1.0,
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

    parser.add_argument('--phase_tumor', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--flanking', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--max_indel_length', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    create_tensor(args)


if __name__ == "__main__":
    main()
