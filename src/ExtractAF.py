import sys
import shlex
import os
import json
import logging
import random
from subprocess import PIPE
from os.path import isfile
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict, OrderedDict

import shared.param as param
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool, vcf_candidates_from
from shared.interval_tree import bed_tree_from, is_region_in
from shared.intervaltree.intervaltree import IntervalTree

logging.basicConfig(format='%(message)s', level=logging.INFO)
BASES = set(list(BASE2NUM.keys()) + ["-"])
no_of_positions = param.no_of_positions
flanking_base_num = param.flankingBaseNum
channel_size = param.channel_size
BASE2NUMBER = dict(zip("ACGTURYSWKMBDHVN-", (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)))
NORMALIZE_NUM = param.NORMALIZE_NUM
MAX_BQ = 40.0
MAX_MQ = 60.0
MAX_AF = 1.0
STRAND_0 = 100
STRAND_1 = 50
HAP_TYPE = dict(zip((1, 0, 2), (30, 60, 90)))  # hap1 UNKNOWN H2
ACGT_NUM = dict(zip("ACGT+-*#N", (100, 25, 75, 50, -50, -100, 0, 0, 100)))


def _normalize_bq(x):
    return int(NORMALIZE_NUM * min(x, MAX_BQ) / MAX_BQ)


def _normalize_mq(x):
    return int(NORMALIZE_NUM * min(x, MAX_MQ) / MAX_MQ)


def _normalize_af(x):
    return int(NORMALIZE_NUM * min(x, MAX_AF) / MAX_AF)


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

    def update_infos(self):
        # only proceed when variant exists in candidate windows which greatly improves efficiency
        self.update_info = True
        self.read_name_dict = dict(zip(self.read_name_list, self.base_list))
        self.mapping_quality = [_normalize_mq(phredscore2raw_score(item)) for item in self.raw_mapping_quality]
        self.base_quality = [_normalize_bq(phredscore2raw_score(item)) for item in self.raw_base_quality]

        for read_name, base_info, bq, mq in zip(self.read_name_list, self.base_list, self.base_quality,
                                                self.mapping_quality):
            read_channel, ins_base, query_base = get_tensor_info(base_info, bq, self.ref_base, mq)
            self.read_info[read_name] = (read_channel, ins_base)


class PhasingRead(object):
    def __init__(self):
        self.read_seq = defaultdict(str)
        self.read_start = None
        self.read_end = None


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


def sorted_by_hap_read_name(center_pos, haplotag_dict, pileup_dict, hap_dict, platform):
    """
    Sort by reads haplotype after haplotag reads otherwise sort by read start position.
    center_pos: define the center candidate position for proccessing.
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
    matrix_depth = param.matrix_depth_dict[platform]
    if len(all_nearby_read_name) > matrix_depth:
        # set same seed for reproducibility
        random.seed(0)
        indices = random.sample(range(len(all_nearby_read_name)), matrix_depth)
        all_nearby_read_name = [all_nearby_read_name[i] for i in sorted(indices)]
    sorted_read_name_list = []
    for order, read_name in enumerate(all_nearby_read_name):
        hap = max(haplotag_dict[read_name], hap_dict[read_name])  # no phasing is 0
        sorted_read_name_list.append((hap, order, read_name))

    sorted_read_name_list = sorted(sorted_read_name_list)
    return sorted_read_name_list


def get_tensor_info(base_info, bq, ref_base, read_mq):
    """
    Create tensor information for each read level position.
    base_info: base information include all alternative bases.
    bq: normalized base quality.
    ref_base: reference_base: upper reference base for cigar calculation.
    read_mq: read mapping quality.
    """

    base, indel = base_info
    ins_base = ""
    query_base = ""
    read_channel = [0] * channel_size
    if base[0] in '*#':
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

    REF_BASE = ACGT_NUM[ref_base]
    if len(indel) and indel[0] in '+-':
        if indel[0] == "+":
            ins_base = indel[1:].upper()
    read_channel[:5] = REF_BASE, ALT_BASE, strand, read_mq, bq
    query_base = "" if base_upper not in "ACGT" else base_upper
    return read_channel, ins_base, query_base


def decode_pileup_bases(pileup_bases, reference_base, minimum_af_for_candidate, has_pileup_candidates, read_name_list, is_tumor):
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
    # if has_pileup_candidates:
    #     return base_list, None, None, None
    pileup_dict = defaultdict(int)
    base_counter = Counter([''.join(item) for item in base_list])
    alt_dict = dict(Counter([''.join(item).upper() for item in base_list]))

    tumor_alt_dict = dict(Counter([''.join(item).upper() for item, read_name in zip(base_list, read_name_list) if read_name.startswith('t')])) if is_tumor else None
    depth = 0
    for read_idx, (key, count) in enumerate(base_counter.items()):
        if key[0].upper() in 'ACGT':
            pileup_dict[key[0].upper()] += count
            depth += count
        if len(key) > 1 and key[1] == '+':
            pileup_dict['I'] += count
        elif len(key) > 1 and key[1] == '-':
            pileup_dict['D'] += count
        # elif len(key) == 1 and key.upper() != reference_base:

    denominator = depth if depth > 0 else 1
    pileup_list = sorted(list(pileup_dict.items()), key=lambda x: x[1], reverse=True)
    af = (float(pileup_list[1][1]) / denominator) if len(pileup_list) > 1 else ((float(pileup_list[0][1]) / denominator) if len(pileup_list) == 1 and pileup_list[0][0] != reference_base else 0.0)
    # pass_af = len(pileup_list) and (pileup_list[0][0] != reference_base or af >= minimum_af_for_candidate)
    pass_af = len(pileup_list) and (af >= minimum_af_for_candidate)
    af = (float(pileup_list[0][1]) / denominator) if len(pileup_list) >= 1 and pileup_list[0][
        0] != reference_base else af

    if not pass_af:
        return base_list, depth, pass_af, af, "", "", ""

    pileup_list = [[item[0], str(round(item[1]/denominator,3))] for item in pileup_list]
    af_infos = ','.join([item[1] for item in pileup_list if item[0] != reference_base])

    alt_list = sorted(list(alt_dict.items()), key=lambda x: x[1], reverse=True)
    alt_list = [[item[0], str(round(item[1]/denominator,3))] for item in alt_list]
    pileup_infos = ' '.join([item[0] + ':' + item[1] for item in alt_list])

    if tumor_alt_dict is not None:
        tumor_alt_list = sorted(list(tumor_alt_dict.items()), key=lambda x: x[1], reverse=True)
        tumor_alt_list = [[item[0], str(round(item[1]/denominator,3))] for item in tumor_alt_list]
        tumor_pileup_infos = ' '.join([item[0] + ':' + item[1] for item in tumor_alt_list])
    else:
        tumor_pileup_infos = ""
    # pileup_list = [[item[0], str(round(item[1]/denominator,3))] for item in pileup_list]
    # af_infos = ','.join([item[1] for item in pileup_list if item[0] != reference_base])
    # pileup_infos = ' '.join([item[0] + ':' + item[1] for item in pileup_list])

    return base_list, depth, pass_af, af, af_infos, pileup_infos, tumor_pileup_infos


def get_alt_info(center_pos, pileup_dict, ref_seq, reference_sequence, reference_start, hap_dict):
    """
    Get alternative information for representation unification, keep all read level alignment information including phasing info.
    center_pos: center position for processing, default window size = no_of_positions = flankingBaseNum + 1 + flankingBaseNum
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    ref_seq: chunked reference sequence in window, start: center pos - flankingBaseNum, end: center + flankingBaseNum + 1.
    reference_sequence: reference sequence index by contig:start-end. 0-based.
    reference_base: upper reference base for cigar calculation.
    reference_start: upper reference base for cigar calculation.
    hap_dict: dictionary (pos: pos info) which keep the read name and haplotype mapping.
    """

    reference_base = ref_seq[flanking_base_num]
    alt_read_name_dict = defaultdict(set)
    depth = 0
    for (base, indel), read_name in zip(pileup_dict[center_pos].base_list, pileup_dict[center_pos].read_name_list):
        if base in "#*":
            alt_read_name_dict['*'].add(read_name)
            depth += 1
            continue
        depth += 1
        if base.upper() == reference_base and indel == '':
            alt_read_name_dict['R'].add(read_name)
        if indel != '':
            if indel[0] == '+':
                indel = 'I' + base.upper() + indel.upper()[1:]
            else:
                del_bases_num = len(indel[1:])
                del_ref_bases = reference_sequence[
                                center_pos - reference_start + 1:center_pos - reference_start + del_bases_num + 1]
                indel = 'D' + del_ref_bases
            alt_read_name_dict[indel].add(read_name)

        if indel == '' and base.upper() != reference_base:
            alt_read_name_dict['X' + base.upper()].add(read_name)

    for alt_type, read_name_set in list(alt_read_name_dict.items()):
        alt_read_name_dict[alt_type] = ' '.join(
            [read_name + '_' + str(hap_dict[read_name]) for read_name in list(read_name_set)])

    alt_info = str(depth) + '\t' + json.dumps(alt_read_name_dict)

    return alt_info


def generate_tensor(ctg_name, center_pos, sorted_read_name_list, pileup_dict, ref_seq, reference_sequence,
                    reference_start, platform, confident_bed_tree, add_no_phasing_data_training, is_tumor,
                    somatic_flag):
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
    if tensor_depth == 0:
        return None, None
    tensor = [[[0] * tensor_shape[2] for _ in range(tensor_shape[1])] for _ in range(tensor_depth)]
    start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num + 1
    insert_tuple = []

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

    # match deletion cases and bed format
    pass_confident_bed = not len(confident_bed_tree) or is_region_in(confident_bed_tree, ctg_name,
                                                                     center_pos - 2,
                                                                     center_pos + max_del_length + 1)
    if not pass_confident_bed:
        return None, None

    for p in range(start_pos, end_pos):
        if p in pileup_dict and not pileup_dict[p].update_info:
            pileup_dict[p].update_infos()
        for read_idx, read_name_info in enumerate(sorted_read_name_list):
            hap, read_order, read_name = read_name_info
            offset = p - start_pos
            if p in pileup_dict and read_name in pileup_dict[p].read_info:
                read_channel, ins_base = pileup_dict[p].read_info[read_name]
                tensor[read_idx][offset] = read_channel
                if ins_base != '' and p < end_pos - 1:
                    insert_tuple.append((read_idx, offset, ins_base, p))

    for read_idx, p, ins_base, cp in insert_tuple:

        for ins_idx in range(min(len(ins_base), no_of_positions - p)):
            tensor[read_idx][ins_idx + p][6] = ACGT_NUM[ins_base[ins_idx]]

    for row_idx, (hap, _, read_name) in enumerate(sorted_read_name_list):
        af_num = 0
        if read_name in pileup_dict[center_pos].read_name_dict:
            base, indel = pileup_dict[center_pos].read_name_dict[read_name]
            base_upper = base.upper()
            if indel != '':
                if indel[0] == '+':
                    insert_str = ('+' + base_upper + indel.upper()[1:])
                    af_num = alt_dict[insert_str] / max(1, float(depth)) if insert_str in alt_dict else af_num
                else:
                    af_num = alt_dict[indel.upper()] / max(1, float(depth)) if indel.upper() in alt_dict else af_num
            elif base.upper() in alt_dict:
                af_num = alt_dict[base_upper] / max(1, float(depth))
        af_num = _normalize_af(af_num) if af_num != 0 else af_num
        # hap_type = HAP_TYPE[hap]
        hap_type = 100 if is_tumor else 50
        for p in range(no_of_positions):
            if tensor[row_idx][p][2] != 0:  # skip all del #*
                tensor[row_idx][p][5] = af_num
                tensor[row_idx][p][7] = hap_type

    alt_info = []
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

    alt_info = str(depth) + '-' + ' '.join([' '.join([item[0], str(item[1])]) for item in alt_info])
    tensor_string_list = [
        " ".join((" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in tensor)]

    # if add_no_phasing_data_training:
    #     all_hap = [item[0] for item in sorted_read_name_list]
    #     # skip if no phased reads exist
    #     if sum(all_hap) != 0:
    #         raw_read_name_index_mapping = [item[1] for item in sorted(
    #             [(item[1], read_idx) for read_idx, item in enumerate(sorted_read_name_list)])]
    #         no_phasing_tensor = [tensor[read_idx] for read_idx in raw_read_name_index_mapping]
    #         for row_idx in range(len(no_phasing_tensor)):
    #             for p in range(no_of_positions):
    #                 if tensor[row_idx][p][7] > 0:
    #                     tensor[row_idx][p][7] = HAP_TYPE[0]
    #
    #         no_phasing_tensor_string = " ".join(
    #             (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in
    #             no_phasing_tensor)
    #         tensor_string_list.append(no_phasing_tensor_string)
    return '\n'.join(["%s\t%d\t%s\t%s\t%s\t%s\t%s" % (
        ctg_name,
        center_pos,
        ref_seq,
        tensor_string,
        alt_info,
        "tumor" if is_tumor else "normal",
        somatic_flag
    ) for tensor_string in tensor_string_list]), alt_info


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def update_hete_ref(pos, reference_sequence, reference_start, extend_bp, alt_base):
    # if need phasing option enables, will store reference squence near hete snp candidate.
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


def CreateTensorFullAlignment(args):
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    full_aln_regions = args.full_aln_regions
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    need_phasing = args.need_phasing
    samtools_execute_command = args.samtools
    output_depth = args.output_depth
    output_alt_info = args.output_alt_info
    bam_file_path = args.bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    is_full_aln_regions_given = full_aln_regions is not None
    phasing_info_in_bam = args.phasing_info_in_bam
    phasing_window_size = args.phasing_window_size
    extend_bp = param.extend_bp
    minimum_af_for_candidate = args.min_af
    minimum_af_for_truth = args.min_truth_af
    platform = args.platform
    store_tumor_infos = args.store_tumor_infos
    alt_fn = args.alt_fn
    extend_bed = args.extend_bed
    is_extend_bed_file_given = extend_bed is not None
    min_mapping_quality = args.minMQ
    min_base_quality = args.minBQ

    vcf_fn = args.vcf_fn
    is_truth_vcf_provided = vcf_fn is not None
    truths_variant_dict = {}
    if is_truth_vcf_provided:
        from shared.vcf import VcfReader
        unified_vcf_reader = VcfReader(vcf_fn=vcf_fn, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        truths_variant_dict = unified_vcf_reader.variant_dict

    global test_pos
    test_pos = 7357230
    if args.test_pos and test_pos:
        platform ="hifi"
        sample_name = 'hg002'
        subsample_ratio = 1000
        ctg_name="chr1"
        ctg = ctg_name[3:]
        minimum_af_for_candidate = 0.08
        minimum_af_for_truth = 0.05
        samtools_execute_command="/autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/samtools"
        bam_file_path="/autofs/bal33/zxzheng/data/bam/pb/hg004/HG004.GRCh38.consensusalignments.bam"
        bam_file_path="/mnt/bal36/zxzheng/somatic/downsample/pair/tumor_chr1_0.1.bam"
        fasta_file_path = "/mnt/bal36/zxzheng/testData/ont/data/GRCh38_no_alt_analysis_set.fasta"
        chunk_num = 100
        chunk_id = None
        extend_bed = None
        phasing_info_in_bam = True
        output_alt_info=True
        output_depth = True
        alt_fn = None
        # unify_repre_fn = "/autofs/bal33/zxzheng/TMP/tmp1"
        ctg_start = test_pos - 1000
        ctg_end = test_pos + 1000
    else:
        test_pos = None

    hete_snp_pos_dict = defaultdict()
    hete_snp_tree = IntervalTree()
    need_phasing_pos_set = set()
    add_read_regions = True
    if full_aln_regions:

        """
        If given full alignment bed regions, all candidate positions will be directly selected from each row, define as 
        'ctg start end', where 0-based center position is the candidate for full alignment calling.
        if 'need_phasing' option enables, full alignment bed regions will also include nearby heterozygous snp candidates for reads
        haplotag, which is faster than whatshap haplotag with more memory occupation.
        """

        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (full_aln_regions)))
        candidate_file_path_output = candidate_file_path_process.stdout

        ctg_start, ctg_end = float('inf'), 0
        for row in candidate_file_path_output:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name: continue
            position = int(row[1]) + 1
            end = int(row[2]) + 1
            ctg_start = min(position, ctg_start)
            ctg_end = max(end, ctg_end)

            if platform == "ilmn":
                continue
            if len(row) > 3:  # hete snp positions
                center_pos = position + extend_bp + 1
                ref_base, alt_base, genotype, phase_set = row[3].split('-')
                hete_snp_pos_dict[center_pos] = Position(pos=center_pos, ref_base=ref_base, alt_base=alt_base,
                                                         genotype=int(genotype), phase_set=phase_set)
                hete_snp_tree.addi(begin=center_pos - extend_bp, end=center_pos + extend_bp + 1)
            else:
                center = position + (end - position) // 2 - 1
                need_phasing_pos_set.add(center)
        candidate_file_path_output.close()
        candidate_file_path_process.wait()

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    if not full_aln_regions and chunk_id is not None:

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


    need_phasing_pos_set = set([item for item in need_phasing_pos_set if item >= ctg_start and item <= ctg_end])
    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    extend_start, extend_end = None, None
    if is_ctg_range_given:
        extend_start = ctg_start - (phasing_window_size if need_phasing else no_of_positions)
        extend_end = ctg_end + (phasing_window_size if need_phasing else no_of_positions)
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
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    # pileup bed first
    read_name_option = ' --output-QNAME' if store_tumor_infos else ' '
    bed_option = ' -l {}'.format(
        extend_bed) if is_extend_bed_file_given and platform != 'ilmn' else ""
    bed_option = ' -l {}'.format(full_aln_regions) if is_full_aln_regions_given and platform != 'ilmn' else bed_option
    flags_option = ' --excl-flags {}'.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    max_depth_option = ' --max-depth {}'.format(args.max_depth) if args.max_depth > 0 else ""
    reads_regions_option = ' -r {}'.format(" ".join(reads_regions)) if add_read_regions else ""
    # print (add_read_regions, ctg_start, ctg_end, reference_start)
    stdin = None if bam_file_path != "PIPE" else sys.stdin
    bam_file_path = bam_file_path if bam_file_path != "PIPE" else "-"
    samtools_command = "{} mpileup  {} --reverse-del".format(samtools_execute_command,
                                                                                        bam_file_path) + \
                       read_name_option + reads_regions_option + phasing_option + mq_option + bq_option + bed_option + flags_option + max_depth_option
    samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command), stdin=stdin)

    if alt_fn:
        output_alt_fn = alt_fn
        alt_fp = open(output_alt_fn, 'w')

    is_tumor = alt_fn.split('/')[-2].startswith('tumor')
    has_pileup_candidates = len(need_phasing_pos_set)
    for row in samtools_mpileup_process.stdout:  # chr position N depth seq BQ read_name mapping_quality phasing_info
        columns = row.strip().split('\t')
        pos = int(columns[1])
        # pos that near bed region should include some indel cover in bed
        # pass_extend_bed = not is_extend_bed_file_given or is_region_in(extend_bed_tree,
        #                                                                ctg_name, pos - 1,
        #                                                                pos + 1)
        # pass_ctg_range = not ctg_start or (pos >= ctg_start and pos <= ctg_end)
        # if not has_pileup_candidates and not pass_extend_bed and pass_ctg_range:
        #     continue
        pileup_bases = columns[4]
        # raw_base_quality = columns[5]
        read_name_list = columns[6].split(',') if store_tumor_infos else []
        # raw_mapping_quality = columns[7]
        reference_base = evc_base_from(reference_sequence[pos - reference_start].upper())  # ev
        is_truth_candidate = pos in truths_variant_dict
        minimum_af = minimum_af_for_truth if is_truth_candidate and minimum_af_for_truth else minimum_af_for_candidate
        base_list, depth, pass_af, af, af_infos, pileup_infos, tumor_pileup_infos = decode_pileup_bases(pileup_bases=pileup_bases,
                                                            reference_base=reference_base,
                                                            minimum_af_for_candidate=minimum_af,
                                                            has_pileup_candidates=has_pileup_candidates,
                                                            read_name_list=read_name_list,
                                                            is_tumor=is_tumor
                                                            )

        # if len(need_phasing_pos_set):
        #     if alt_fn and pos in need_phasing_pos_set:
        #         alt_fp.write('\t'.join([ctg_name + '\t' + str(pos), str(af)]) + '\n')
        if pass_af and alt_fn:
            depth_list = [str(depth)] if output_depth else []
            alt_info_list = [af_infos, pileup_infos, tumor_pileup_infos] if output_alt_info else []
            alt_fp.write('\t'.join([ctg_name, str(pos), reference_base] + depth_list + alt_info_list) + '\n')
    samtools_mpileup_process.stdout.close()
    samtools_mpileup_process.wait()

    if alt_fn:
        alt_fp.close()


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.00,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--min_truth_af', type=float, default=None,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for advanced users
    parser.add_argument('--minCoverage', type=float, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$minMQ are filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="EXPERIMENTAL: If set, bases with base quality with <$minBQ are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=144,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--output_depth', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--output_alt_info', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--phasing_window_size', type=int, default=param.phasing_window_size,
                        help="DEBUG: The window size for read phasing")

    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--alt_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")

    parser.add_argument('--store_tumor_infos', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")

    # options for internal process control
    ## Path to the 'zstd' compression
    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help=SUPPRESS)

    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=str2bool, default=1,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Only call variant in phased vcf file
    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Apply no phased data in training. Only works in data training, default: False
    parser.add_argument('--add_no_phasing_data_training', action='store_true',
                        help=SUPPRESS)

    ## Output representation unification infos, which refines training labels
    parser.add_argument('--unify_repre', action='store_true',
                        help=SUPPRESS)

    ## Path of representation unification output
    parser.add_argument('--unify_repre_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help=SUPPRESS)

    ## Use Clair3's own phasing module for read level phasing when creating tensor, compared to using Whatshap, speed is faster but has higher memory footprint, default: False
    parser.add_argument('--need_phasing', action='store_true',
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform. Greatly boost indel performance in trade of running time
    parser.add_argument('--need_realignment', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    CreateTensorFullAlignment(args)


if __name__ == "__main__":
    main()
