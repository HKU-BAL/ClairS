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
from subprocess import PIPE

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
base_index = dict(zip(phase_channel, list(range(16))))

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


def decode_pileup_bases(pos, pileup_bases, reference_base,  minimum_snp_af_for_candidate, minimum_indel_af_for_candidate, \
                        has_pileup_candidates, candidates_type_dict,is_tumor, mapping_quality, base_quality, base_list=None, \
                        phasing_info=None, platform="ont"):
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

    pileup_tensor = [0] * (channel_size if phasing_info is None else (channel_size + len(phase_channel)))
    is_candidate = pos in candidates_type_dict
    if base_list is None:
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
    alt_info_dict = defaultdict(int)
    for key, count in base_counter.items():
        if len(key) == 1:
            if key.upper() in 'ACGT':
                pileup_dict[key.upper()] += count
                if is_candidate and key.upper() != reference_base:
                    alt_info_dict['X' + key.upper()] += count
                depth += count
                pileup_tensor[BASE2INDEX[key]] += count
            elif key in '#*':
                pileup_tensor[BASE2INDEX[key]] += count
                depth += count
        elif key[1] == '+':
            pileup_dict['I'] += count
            if is_candidate:
                alt_info_dict['I' + key[0].upper() + key[2:].upper()] += count
            # two strand
            if key[0] in 'ACGTN*':
                pileup_tensor[BASE2INDEX["I"]] += count
                max_ins_0 = max(max_ins_0, count)
            else:
                pileup_tensor[BASE2INDEX["i"]] += count
                max_ins_1 = max(max_ins_1, count)
        elif key[1] == '-':
            pileup_dict['D'] += count
            if is_candidate:
                alt_info_dict['D' + key[0].upper() + key[2:].upper()] += count
            max_del_length = max(max_del_length, len(key[1:]))
            # two strand
            if key[0] in 'N*ACGT':
                pileup_tensor[BASE2INDEX["D"]] += count
                max_del_0 = max(max_del_0, count)
            else:
                pileup_tensor[BASE2INDEX["d"]] += count
                max_del_1 = max(max_del_1, count)
    alt_info = str(depth) + '-' + ' '.join(
                [' '.join([item[0], str(item[1])]) for item in alt_info_dict.items()]) + '-'
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


def find_tumor_alt_match(center_pos,
                         sorted_read_name_list,
                         read_name_dict,
                         all_nearby_read_name,
                         truths_variant_dict,
                         normal_output_bam_prefix= 'n',
                         tumor_output_bam_prefix='t',
                         proportion=None):
    # if proportion is not None and float(proportion) == 1.0:
    #     # all reads are from tumor reads
    #     tumor_reads = all_nearby_read_name
    #     normal_reads = []
    # else:
    tumor_reads = [read_name for read_name in all_nearby_read_name if read_name.startswith(tumor_output_bam_prefix)]
    normal_reads = [read_name for read_name in all_nearby_read_name if read_name.startswith(normal_output_bam_prefix)]
    ref_base, alt_base = truths_variant_dict[center_pos].reference_bases, truths_variant_dict[center_pos].alternate_bases[0]
    is_ins = len(alt_base) > 1 and len(ref_base) == 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    is_snp = len(ref_base) == 1 and len(alt_base) == 1

    matched_read_name_set = set()
    normal_read_name_set = set(normal_reads)
    tumor_read_name_set = set(tumor_reads)
    for read_name in tumor_reads:
        if read_name in read_name_dict:
            base, indel = read_name_dict[read_name]
            base_upper = base.upper()
            if is_ins and base_upper + indel[1:].upper() == alt_base:
                matched_read_name_set.add(read_name)
            elif is_del and len(indel[1:].upper()) == len(ref_base[1:]):
                matched_read_name_set.add(read_name)
            elif is_snp and base_upper == alt_base:
                matched_read_name_set.add(read_name)
    return matched_read_name_set, normal_read_name_set, tumor_read_name_set


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
    extend_bp = param.extend_bp

    minimum_snp_af_for_candidate = args.snv_min_af
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
    phasing_info_in_bam = args.phase_tumor and tensor_sample_mode and args.platform != 'ilmn'
    global test_pos
    test_pos = None
    candidates_pos_set = set()
    candidates_type_dict = defaultdict(str)
    add_read_regions = True
    training_mode = args.add_phasing_info or args.training_mode
    truth_vcf_fn = args.truth_vcf_fn
    is_truth_vcf_provided = truth_vcf_fn is not None
    truths_variant_dict = {}
    proportion = args.proportion

    normal_output_bam_prefix = args.normal_output_bam_prefix
    tumor_output_bam_prefix = args.tumor_output_bam_prefix

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
            # if variant_type == 'hem'
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

        if bam_file_path != "PIPE":
            bam_file_path += '.{}_{}'.format(ctg_start, ctg_end)
            add_read_regions = False
        if bam_file_path == "PIPE":
            add_read_regions = False

    # preparation for candidates near variants
    candidates_pos_set = set([item for item in candidates_pos_set if item >= ctg_start and item <= ctg_end])
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
    samtools_mpileup_min_mq = 0
    mq_option = ' --min-MQ {}'.format(samtools_mpileup_min_mq)
    output_mq, output_read_name = True, training_mode
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
    stdin = None if bam_file_path != "PIPE" else sys.stdin
    bam_file_path = bam_file_path if bam_file_path != "PIPE" else "-"
    samtools_command = "{} mpileup {} --reverse-del".format(samtools_execute_command, bam_file_path) + \
                       output_read_name_option + output_mq_option + reads_regions_option + phasing_option + mq_option + bq_option + bed_option + flags_option + max_depth_option
    samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command), stdin=stdin, stderr=PIPE)

    is_tumor = "tumor_" in bam_file_path or tensor_sample_mode
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
    extend_bp_distance = no_of_positions + param.extend_bp


    extend_bed_tree = bed_tree_from(bed_file_path=extend_bed,
                                    contig_name=ctg_name,
                                    bed_ctg_start=extend_start,
                                    bed_ctg_end=extend_end)

    pileup_tensors = [[0] * channel_size] * (extend_end - extend_start)

    alt_info_dict = defaultdict()
    pileup_info_dict = defaultdict()
    candidate_pos_list = sorted(list(candidates_pos_set))
    has_pileup_candidates = len(candidates_pos_set)
    use_alt_base = False if platform != 'ilmn' else param.use_alt_base
    hap_dict = defaultdict(int)

    def samtools_pileup_generator_from(samtools_mpileup_process):
        current_pos_index = 0
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
            read_name_list = columns[7].split(',') if training_mode else []
            reference_base = reference_sequence[pos - reference_start].upper()
            if reference_base not in 'ACGT':
                continue

            mapping_quality = [ord(mq) - 33 for mq in raw_mapping_quality]
            base_quality = [ord(mq) - 33 for mq in raw_base_quality]

            if phasing_info_in_bam:
                phasing_info = columns[8].split(',')
                for hap_idx, hap in enumerate(phasing_info):
                    if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                        hap_dict[read_name_list[hap_idx]] = int(hap)

            else:
                phasing_info = None
            pileup_tensor, base_list, depth, pass_af, af, alt_info = decode_pileup_bases(pos=pos,
                                                                pileup_bases=pileup_bases,
                                                                reference_base=reference_base,
                                                                minimum_snp_af_for_candidate=minimum_snp_af_for_candidate,
                                                                minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
                                                                has_pileup_candidates=has_pileup_candidates,
                                                                candidates_type_dict=candidates_type_dict,
                                                                mapping_quality=mapping_quality,
                                                                base_quality=base_quality,
                                                                phasing_info=phasing_info,
                                                                # use_tensor_sample_mode=use_tensor_sample_mode,
                                                                is_tumor=is_tumor)
            if training_mode:
                for b_idx, base in enumerate(base_list):
                    if base[0] == '#' or (base[0] >= 'a' and base[0] <= 'z'):
                        read_name_list[b_idx] += '_1'  # reverse
                    else:
                        read_name_list[b_idx] += '_0'  # forward

                pileup_info_dict[pos] = (base_list, read_name_list, mapping_quality, base_quality, phasing_info)
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
                yield candidate_pos_list[current_pos_index]
                current_pos_index += 1

        while current_pos_index != len(candidate_pos_list):
            yield candidate_pos_list[current_pos_index]

            current_pos_index += 1

    samtools_pileup_generator = samtools_pileup_generator_from(samtools_mpileup_process)

    for pos in samtools_pileup_generator:
        ref_seq = reference_sequence[
                  pos - reference_start - flanking_base_num: pos - reference_start + flanking_base_num + 1].upper()

        start_index = pos - flanking_base_num - extend_start
        end_index = start_index + no_of_positions
        if start_index < 0 or end_index >= extend_end - extend_start:
            continue
        if sum([1 for tensor in pileup_tensors[start_index:end_index] if len(tensor) == 0]) > 0:
            continue

        use_tensor_sample_mode = tensor_sample_mode and (
                candidates_type_dict[pos] == 'homo_somatic' or candidates_type_dict[
            pos] == 'hetero_somatic') and pos in truths_variant_dict and training_mode

        tensor_string_list = []
        alt_info_list = []
        if pos not in alt_info_dict:
            continue
        if use_tensor_sample_mode:

            base_list, read_name_list = pileup_info_dict[pos][0], pileup_info_dict[pos][1]
            read_name_dict = dict(zip(read_name_list, base_list))
            center_pos = pos
            all_nearby_read_name = []
            start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num + 1
            for p in range(start_pos, end_pos):

                if p in pileup_info_dict:
                    all_nearby_read_name += pileup_info_dict[p][1]
            all_nearby_read_name = list(OrderedDict.fromkeys(all_nearby_read_name))  # have sorted by order
            tumor_reads_meet_alt_info_set, normal_read_name_set, tumor_read_name_set = find_tumor_alt_match(center_pos,
                                                                                                            read_name_list,
                                                                                                            read_name_dict,
                                                                                                            all_nearby_read_name,
                                                                                                            truths_variant_dict,
                                                                                                            normal_output_bam_prefix=normal_output_bam_prefix,
                                                                                                            tumor_output_bam_prefix=tumor_output_bam_prefix,
                                                                                                            proportion=proportion)

            if len(tumor_reads_meet_alt_info_set) == 0:
                print("No reads support tumor alternative in pos:{}".format(center_pos))
                continue
            min_tumor_support_read_num = param.min_tumor_support_read_num
            tumor_read_name_list = list(tumor_read_name_set)
            normal_read_name_list = list(normal_read_name_set)
            tumor_reads_num = len(tumor_read_name_set)
            tumor_reads_meet_alt_info_num = len(tumor_reads_meet_alt_info_set)
            tumor_read_porportion = tumor_reads_meet_alt_info_num / float(tumor_reads_num)
            sampled_reads_num_list = []

            partition = 3
            if param.use_beta_subsampling:
                beta_acc_per = param.beta_acc_per
                sampled_reads_num_list.append(len(tumor_read_name_list))
                for s_idx in range(partition):
                    random.seed(pos + s_idx)
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

            elif param.use_exp_subsampling:
                for p in range(1, partition):
                    af = 1 / float(2 ** p)
                    sampled_read_num = int((len(tumor_read_name_set) * af) / tumor_read_porportion)
                    if sampled_read_num >= min_tumor_support_read_num and sampled_read_num <= len(tumor_read_name_set):
                        sampled_reads_num_list.append(sampled_read_num)


            for read_num in sampled_reads_num_list:
                random.seed(len(tumor_read_name_list) + read_num)
                sampled_tumor_read_name_list = random.sample(tumor_read_name_list, read_num)
                sampled_tumor_read_name_meet_alt_set = set(sampled_tumor_read_name_list).intersection(
                    tumor_reads_meet_alt_info_set)
                # in training mode, the candidates should have enough read support
                if len(sampled_tumor_read_name_meet_alt_set) < min_tumor_support_read_num:
                    continue

                if len(sampled_tumor_read_name_list) + len(normal_read_name_list) < len(all_nearby_read_name) * 0.3 and not use_alt_base:
                    continue

                if phasing_info_in_bam:
                    hap_list = [0, 0, 0]
                    for rn in sampled_tumor_read_name_meet_alt_set:
                        if rn[:-2] in hap_dict:
                            hap_list[hap_dict[rn[:-2]]] += 1
                    if hap_list[1] > 0 and hap_list[2] > 0:
                        print("pos exist in both side")
                        continue

                # check the strand info of sampled_tumor_read_name_list
                forward_strand = sum([1 for rn in sampled_tumor_read_name_meet_alt_set if rn.endswith('1')])
                if forward_strand == 0 or forward_strand == len(sampled_tumor_read_name_meet_alt_set):
                    continue

                if use_alt_base:
                    sorted_read_name_set = set(normal_read_name_list + tumor_read_name_list)
                else:
                    sorted_read_name_set = set(normal_read_name_list + sampled_tumor_read_name_list)
                tmp_pileup_tensors = []
                tmp_alt_info = []
                for p in range(start_pos, end_pos):
                    pre_base_list, pre_read_name_list, pre_mapping_quality, pre_base_quality, pre_phasing_info = pileup_info_dict[p]
                    base_list, read_name_list, mapping_quality, base_quality = [], [], [], []
                    phasing_info = [] if phasing_info_in_bam else None
                    for b_idx, (base, read_name, mq, bq) in enumerate(zip(pre_base_list, pre_read_name_list, pre_mapping_quality, pre_base_quality)):
                        if read_name in sorted_read_name_set:
                            if use_alt_base and p == pos and read_name not in sampled_tumor_read_name_meet_alt_set:
                                if base[0].upper() == truths_variant_dict[pos].alternate_bases[0]:
                                    base = [ref_seq[flanking_base_num] if base[0] in "ACGT" else ref_seq[flanking_base_num].lower(), ""]
                            base_list.append(base)
                            read_name_list.append(read_name)
                            mapping_quality.append(mq)
                            base_quality.append(bq)
                            if phasing_info_in_bam:
                                phasing_info.append(pre_phasing_info[b_idx])

                    reference_base = reference_sequence[p - reference_start].upper()
                    pileup_tensor, base_list, depth, pass_af, af, alt_info = decode_pileup_bases(pos=p,
                                                                                                 pileup_bases="",
                                                                                                 reference_base=reference_base,
                                                                                                 minimum_snp_af_for_candidate=minimum_snp_af_for_candidate,
                                                                                                 minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
                                                                                                 has_pileup_candidates=has_pileup_candidates,
                                                                                                 candidates_type_dict=candidates_type_dict,
                                                                                                 mapping_quality=mapping_quality,
                                                                                                 base_quality=base_quality,
                                                                                                 base_list=base_list,
                                                                                                 phasing_info=phasing_info,
                                                                                                 is_tumor=is_tumor)
                    tmp_pileup_tensors.append(pileup_tensor)
                    if p == pos:
                        tmp_alt_info = alt_info
                if sum([1 for tensor in tmp_pileup_tensors if len(tensor) == 0]) > 0:
                    continue

                tensor_string_list.append(" ".join(" ".join("%d" % x for x in innerlist) for innerlist in tmp_pileup_tensors))
                alt_info_list.append(tmp_alt_info)
                if 0:
                    import numpy as np
                    a = np.array(tmp_pileup_tensors)
                    b = np.array(pileup_tensors[start_index:end_index])
        else:
            tensor_string_list.append(" ".join(" ".join("%d" % x for x in innerlist) for innerlist in pileup_tensors[start_index:end_index]))
            alt_info_list.append(alt_info_dict[pos])
        if len(tensor_string_list) > 0:
            variant_type = candidates_type_dict[pos] if pos in candidates_type_dict else 'unknown'
            tensor = '\n'.join(["%s\t%d\t%s\t%s\t%s\t%s\t%s" % (
                ctg_name,
                pos,
                ref_seq,
                tensor_string,
                alt_info,
                "tumor" if is_tumor else "normal",
                variant_type
            ) for tensor_string, alt_info in zip(tensor_string_list, alt_info_list)])
            tensor_can_fp.stdin.write(tensor+"\n")
            if 0:
                import numpy as np
                a = []
                for tensor_str in tensor_string_list:
                    a += tensor_str.split(' ')
                a = np.array(a, dtype=float).reshape((-1, no_of_positions, channel_size + 16 if phasing_info_in_bam else channel_size))

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
    parser = ArgumentParser(description="Create pileup tesnors for training")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
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
    parser.add_argument('--indel_min_af', type=float, default=0.2,
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

    parser.add_argument('--tensor_sample_mode', type=str2bool, default=1,
                        help=SUPPRESS)

    parser.add_argument('--phase_tumor', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--training_mode', type=str2bool, default=0,
                        help=SUPPRESS)

    parser.add_argument('--add_phasing_info', type=str2bool, default=0,
                        help=SUPPRESS)
    parser.add_argument('--proportion', type=float, default=1.0,
                        help=SUPPRESS)

    parser.add_argument('--normal_output_bam_prefix', type=str, default='n',
                        help="Normal output BAM prefix")

    parser.add_argument('--tumor_output_bam_prefix', type=str, default='t',
                        help="Tumor output BAM prefix")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    create_tensor(args)


if __name__ == "__main__":
    main()
