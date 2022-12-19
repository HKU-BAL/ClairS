import sys
import shlex
import os
import logging
import subprocess

from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter
from shared.utils import subprocess_popen, file_path_from, region_from, \
    reference_sequence_from, str2bool
from shared.interval_tree import bed_tree_from, is_region_in

logging.basicConfig(format='%(message)s', level=logging.INFO)

def decode_pileup_bases(pileup_bases, reference_base, min_coverage, minimum_snv_af_for_candidate,
                        minimum_indel_af_for_candidate, alternative_base_num, has_pileup_candidates, read_name_list,
                        is_tumor, platform="ont"):
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

    pileup_dict = defaultdict(int)
    base_counter = Counter([''.join(item) for item in base_list])
    alt_dict = dict(Counter([''.join(item).upper() for item in base_list]))

    tumor_alt_dict = dict(Counter([''.join(item).upper() for item, read_name in zip(base_list, read_name_list) if
                                   read_name.startswith('t')])) if is_tumor else None
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

    denominator = depth if depth > 0 else 1
    pileup_list = sorted(list(pileup_dict.items()), key=lambda x: x[1], reverse=True)

    pass_snv_af = False
    pass_indel_af = False

    pass_depth = depth > min_coverage
    for item, count in pileup_list:
        if item == reference_base:
            continue
        elif item[0] in 'ID':
            pass_indel_af = (pass_indel_af or (float(count) / denominator >= minimum_indel_af_for_candidate))
            continue
        pass_snv_af = pass_snv_af or (float(count) / denominator >= minimum_snv_af_for_candidate) and (
                    alternative_base_num is not None and count >= alternative_base_num)


    af = (float(pileup_list[1][1]) / denominator) if len(pileup_list) > 1 else 0.0
    af = (float(pileup_list[0][1]) / denominator) if len(pileup_list) >= 1 and pileup_list[0][
        0] != reference_base else af

    pass_af = (pass_snv_af or pass_indel_af) and pass_depth

    alt_list = sorted(list(alt_dict.items()), key=lambda x: x[1], reverse=True)
    alt_list = [[item[0], str(round(item[1] / denominator, 3))] for item in alt_list if
                item[0].upper() != reference_base]

    if not pass_af:
        return base_list, depth, pass_af, af, "", "", "", alt_list

    pileup_list = [[item[0], str(round(item[1] / denominator, 3))] for item in pileup_list]
    af_infos = ','.join([item[1] for item in pileup_list if item[0] != reference_base])

    pileup_infos = ' '.join([item[0] + ':' + item[1] for item in alt_list])

    if tumor_alt_dict is not None:
        tumor_alt_list = sorted(list(tumor_alt_dict.items()), key=lambda x: x[1], reverse=True)
        tumor_alt_list = [[item[0], str(round(item[1] / denominator, 3))] for item in tumor_alt_list]
        tumor_pileup_infos = ' '.join([item[0] + ':' + item[1] for item in tumor_alt_list])
    else:
        tumor_pileup_infos = ""

    return base_list, depth, pass_af, af, af_infos, pileup_infos, tumor_pileup_infos, alt_list


def extract_pair_candidates(args):
    ctg_start = args.ctg_start
    ctg_end = args.ctg_end
    fasta_file_path = args.ref_fn
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    output_depth = args.output_depth
    output_alt_info = args.output_alt_info
    tumor_bam_file_path = args.tumor_bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    minimum_snv_af_for_candidate = args.snv_min_af
    normal_snv_max_af = args.normal_snv_max_af
    minimum_indel_af_for_candidate = args.indel_min_af
    minimum_snv_af_for_truth = args.min_truth_snv_af
    minimum_indel_af_for_truth = args.min_truth_snv_af
    alternative_base_num = args.alternative_base_num
    split_bed_size = param.split_bed_size
    candidates_folder = args.candidates_folder
    min_coverage = args.min_coverage
    platform = args.platform
    store_tumor_infos = args.store_tumor_infos
    alt_fn = args.alt_fn
    confident_bed_fn = file_path_from(args.bed_fn, allow_none=True, exit_on_not_found=False)
    is_confident_bed_file_given = confident_bed_fn is not None
    min_mapping_quality = args.min_mq
    min_base_quality = args.min_bq if args.min_bq is not None else param.min_bq_dict[platform]
    flankingBaseNum = param.flankingBaseNum if args.flanking is None else args.flanking
    no_of_positions = 2 * flankingBaseNum + 1
    vcf_fn = args.vcf_fn
    genotyping_mode = vcf_fn is not None
    truth_vcf_fn = args.truth_vcf_fn
    is_truth_vcf_provided = truth_vcf_fn is not None

    truths_variant_dict = {}
    if is_truth_vcf_provided:
        unified_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        truths_variant_dict = unified_vcf_reader.variant_dict

    candidates_pos_set = set()
    add_read_regions = True

    if genotyping_mode:
        vcf_reader = VcfReader(vcf_fn=vcf_fn, ctg_name=ctg_name, is_var_format=False)
        vcf_reader.read_vcf()
        variant_dict = vcf_reader.variant_dict

        candidates_list = sorted(list(variant_dict.keys()))
        if candidates_folder is not None and len(candidates_list):
            candidates_regions = []
            region_num = len(candidates_list) // split_bed_size + 1 if len(
                candidates_list) % split_bed_size else len(candidates_list) // split_bed_size

            for idx in range(region_num):
                # a windows region for create tensor # samtools mpileup not include last position
                split_output = candidates_list[idx * split_bed_size: (idx + 1) * split_bed_size]
                output_path = os.path.join(candidates_folder, '{}.{}_{}_{}'.format(ctg_name, chunk_id, idx, region_num))
                candidates_regions.append(output_path)
                with open(output_path, 'w') as output_file:
                    output_file.write('\n'.join(
                        ['\t'.join([ctg_name, str(x - flankingBaseNum - 1), str(x + flankingBaseNum + 1)]) for x in
                         split_output]) + '\n')  # bed format

            all_candidates_regions_path = os.path.join(candidates_folder,
                                                     'CANDIDATES_FILE_{}_{}'.format(ctg_name, chunk_id))
            with open(all_candidates_regions_path, 'w') as output_file:
                output_file.write('\n'.join(candidates_regions) + '\n')
        return

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    if chunk_id is not None:

        """
        Whole genome calling option, acquire contig start end position from reference fasta index(.fai), then split the
        reference accroding to chunk id and total chunk numbers.
        """
        if is_confident_bed_file_given:
            # consistent with pileup generation, faster to extract tensor using bed region
            tree, bed_start, bed_end = bed_tree_from(bed_file_path=confident_bed_fn,
                                                     contig_name=ctg_name,
                                                     return_bed_region=True)

            chunk_size = (bed_end - bed_start) // chunk_num + 1 if (bed_end - bed_start) % chunk_num else (
                                                                                                                      bed_end - bed_start) // chunk_num
            ctg_start = bed_start + 1 + chunk_size * chunk_id  # 0-base to 1-base
            ctg_end = ctg_start + chunk_size
        else:
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

    candidates_pos_set = set([item for item in candidates_pos_set if item >= ctg_start and item <= ctg_end])
    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    if is_ctg_range_given:
        extend_start = max(ctg_start - ( no_of_positions), 1)
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

    mq_option = ' --min-MQ {}'.format(min_mapping_quality)
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    read_name_option = ' --output-QNAME' if store_tumor_infos else ' '
    bed_option = ' -l {}'.format(
        confident_bed_fn) if is_confident_bed_file_given else ""
    flags_option = ' --excl-flags {} '.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    # add max-depth cut-off? TODO
    max_depth_option = ' --max-depth {} '.format(args.max_depth) if args.max_depth is not None else " "
    reads_regions_option = ' -r {}'.format(" ".join(reads_regions)) if add_read_regions else ""
    # print (add_read_regions, ctg_start, ctg_end, reference_start)
    stdin = None if tumor_bam_file_path != "PIPE" else sys.stdin
    tumor_bam_file_path = tumor_bam_file_path if tumor_bam_file_path != "PIPE" else "-"
    samtools_command = samtools_execute_command + " mpileup --reverse-del" + read_name_option + reads_regions_option + \
                       mq_option + bq_option + bed_option + flags_option + max_depth_option

    samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command + ' ' +tumor_bam_file_path), stdin=stdin)

    if alt_fn:
        output_alt_fn = alt_fn
        alt_fp = open(output_alt_fn, 'w')

    is_tumor = alt_fn.split('/')[-2].startswith('tumor') if alt_fn else False
    has_pileup_candidates = len(candidates_pos_set)
    candidates_list = []

    candidates_dict = defaultdict(str)
    for row in samtools_mpileup_process.stdout:  # chr position N depth seq BQ read_name mapping_quality phasing_info
        columns = row.strip().split('\t')
        pos = int(columns[1])

        pileup_bases = columns[4]
        read_name_list = columns[6].split(',') if store_tumor_infos else []
        reference_base = reference_sequence[pos - reference_start].upper()
        if reference_base.upper() not in "ACGT":
            continue
        is_truth_candidate = pos in truths_variant_dict
        minimum_snv_af_for_candidate = minimum_snv_af_for_truth if is_truth_candidate and minimum_snv_af_for_truth else minimum_snv_af_for_candidate
        minimum_indel_af_for_candidate = minimum_indel_af_for_truth if is_truth_candidate and minimum_indel_af_for_truth else minimum_indel_af_for_candidate
        base_list, depth, pass_af, af, af_infos, pileup_infos, tumor_pileup_infos, alt_list = decode_pileup_bases(
            pileup_bases=pileup_bases,
            reference_base=reference_base,
            min_coverage=min_coverage,
            minimum_snv_af_for_candidate=minimum_snv_af_for_candidate,
            minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
            alternative_base_num=alternative_base_num,
            has_pileup_candidates=has_pileup_candidates,
            read_name_list=read_name_list,
            is_tumor=is_tumor
        )

        if pass_af and alt_fn:
            depth_list = [str(depth)] if output_depth else []
            alt_info_list = [af_infos, pileup_infos, tumor_pileup_infos] if output_alt_info else []
            alt_fp.write('\t'.join([ctg_name, str(pos), reference_base] + depth_list + alt_info_list) + '\n')

        if pass_af:
            candidates_list.append(pos)
            candidates_dict[pos] = (alt_list, depth)

    # scan the normal_bam
    bed_path = os.path.join(candidates_folder, "bed", '{}_{}.bed'.format(ctg_name, chunk_id))
    if not os.path.exists(os.path.join(candidates_folder, 'bed')):
        output = subprocess.run("mkdir -p {}".format(os.path.join(candidates_folder, 'bed')), shell=True)
    output_bed = open(bed_path, 'w')
    for pos in sorted(candidates_list):
        output_bed.write('\t'.join([ctg_name, str(pos - 1), str(pos)]) + '\n')
    output_bed.close()

    candidates_set = set(candidates_list)
    normal_samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command + ' ' +args.normal_bam_fn + ' -l ' + bed_path), stdin=stdin)

    high_normal_af_set = set()
    high_af_gap_set = set()
    for row in normal_samtools_mpileup_process.stdout:  # chr position N depth seq BQ read_name mapping_quality phasing_info
        columns = row.strip().split('\t')
        pos = int(columns[1])

        pileup_bases = columns[4]
        read_name_list = columns[6].split(',') if store_tumor_infos else []
        reference_base = reference_sequence[pos - reference_start].upper()
        if pos not in candidates_set:
            continue
        if reference_base.upper() not in "ACGT":
            continue
        is_truth_candidate = pos in truths_variant_dict
        minimum_snv_af_for_candidate = minimum_snv_af_for_truth if is_truth_candidate and minimum_snv_af_for_truth else minimum_snv_af_for_candidate
        minimum_indel_af_for_candidate = minimum_indel_af_for_truth if is_truth_candidate and minimum_indel_af_for_truth else minimum_indel_af_for_candidate
        base_list, depth, pass_af, af, af_infos, pileup_infos, normal_pileup_infos, normal_alt_list = decode_pileup_bases(
            pileup_bases=pileup_bases,
            reference_base=reference_base,
            min_coverage=min_coverage,
            minimum_snv_af_for_candidate=minimum_snv_af_for_candidate,
            minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
            alternative_base_num=alternative_base_num,
            has_pileup_candidates=has_pileup_candidates,
            read_name_list=read_name_list,
            is_tumor=is_tumor
        )

        tumor_alt_list, tumor_depth = candidates_dict[pos]
        tumor_info = [item for item in tumor_alt_list if item[0] in "ACGT"]
        if len(tumor_info) == 0:
            candidates_set.remove(pos)
            continue
        alt_base, tumor_af = tumor_info[0]
        normal_info = [item for item in normal_alt_list if item[0] == alt_base]
        tumor_af = float(tumor_af)
        if depth <= min_coverage:
            continue
        if len(normal_info) == 0:
            continue
        normal_af = float(normal_info[0][1])
        if normal_af > 0 and tumor_af > 0:
            if normal_af >= normal_snv_max_af:
                if not (tumor_af >= normal_af * 6):
                    candidates_set.remove(pos)
                    high_normal_af_set.add(pos)
            elif tumor_af <= normal_af * 4 and tumor_af != 0:
                candidates_set.remove(pos)
                high_af_gap_set.add(pos)

    candidates_list = [pos for pos in candidates_list if pos in candidates_set]

    gen_vcf = False
    if gen_vcf:
        truth_vcf_fn = args.truth_vcf_fn
        truth_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, show_ref=False, keep_row_str=True,
                                     skip_genotype=True)
        truth_vcf_reader.read_vcf()
        truth_variant_dict = truth_vcf_reader.variant_dict

        vcf_writer = VcfWriter(vcf_fn=os.path.join(candidates_folder, "{}_{}.vcf".format(ctg_name, chunk_id)),
                               ref_fn=fasta_file_path, ctg_name=ctg_name, show_ref_calls=True)
        for pos in candidates_list:
            genotype = '1/1'
            ref_base, alt_base = "A", "A"
            if pos in truth_variant_dict:
                print(ctg_name, pos, "in truth set")
                continue

            tumor_alt_list, tumor_depth = candidates_dict[pos]
            tumor_info = [item for item in tumor_alt_list if item[0] in "ACGT"]
            if len(tumor_info) == 0:
                candidates_set.remove(pos)
                print(pos, "gen vcf not found tumor")
                continue
            alt_base, tumor_af = tumor_info[0]
            ref_base = reference_sequence[pos - reference_start].upper()
            vcf_writer.write_row(POS=pos,
                                 REF=ref_base,
                                 ALT=alt_base,
                                 QUAL=10,
                                 GT=genotype,
                                 DP=10,
                                 AF=float(tumor_af))
        vcf_writer.close()

    print("[INFO] {} high_normal_af_count/high_af_gap_set: {}/{}".format(ctg_name, len(high_normal_af_set),
                                                                         len(high_af_gap_set)))

    if candidates_folder is not None and len(candidates_list):
        all_candidates_regions = []
        region_num = len(candidates_list) // split_bed_size + 1 if len(
            candidates_list) % split_bed_size else len(candidates_list) // split_bed_size

        for idx in range(region_num):
            # a windows region for create tensor # samtools mpileup not include last position
            split_output = candidates_list[idx * split_bed_size: (idx + 1) * split_bed_size]
            output_path = os.path.join(candidates_folder, '{}.{}_{}_{}'.format(ctg_name, chunk_id, idx, region_num))
            all_candidates_regions.append(output_path)
            with open(output_path, 'w') as output_file:
                output_file.write('\n'.join(
                    ['\t'.join([ctg_name, str(x - flankingBaseNum - 1), str(x + flankingBaseNum + 1)]) for x in
                     split_output]) + '\n')  # bed format

        all_candidates_regions_path = os.path.join(candidates_folder,
                                                   'CANDIDATES_FILE_{}_{}'.format(ctg_name, chunk_id))
        with open(all_candidates_regions_path, 'w') as output_file:
            output_file.write('\n'.join(all_candidates_regions) + '\n')

    samtools_mpileup_process.stdout.close()
    samtools_mpileup_process.wait()

    if alt_fn:
        alt_fp.close()


def main():
    parser = ArgumentParser(description="Generate normal-tumor pair variant candidate for tensor creation")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input, default: %(default)s")

    parser.add_argument('--candidates_folder', type=str, default=None,
                        help="Output candidate folder to store the candidate bed information, required")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file, default: %(default)s")

    parser.add_argument('--snv_min_af', type=float, default=param.snv_min_af,
                        help="Minimum SNV allele frequency for a site to be considered as a candidate site in tumor sample, default: %(default)f")

    parser.add_argument('--normal_snv_max_af', type=float, default=param.normal_snv_max_af,
                        help="Maximum SNV allele frequency for a site to be considered as a somatic candidate site in normal sample, default: %(default)f")

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
                        help="EXPERIMENTAL: Maximum depth to be processed. default: %(default)s")

    parser.add_argument('--alternative_base_num', type=int, default=param.alternative_base_num,
                        help="EXPERIMENTAL: Minimum alternative base number to process a candidate. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--output_depth', type=str2bool, default=False,
                        help="DEBUG: Output the candidate coverage for debugging, default: False")

    parser.add_argument('--output_alt_info', type=str2bool, default=False,
                        help="DEBUG: Output the candidate alternative information for debugging, default: False")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation")

    parser.add_argument('--max_pro', type=int, default=None,
                        help="DEBUG: The window size for read phasing")

    # options for internal process control
    ## Minimum SNV allele frequency for a truth site for training
    parser.add_argument('--min_truth_snv_af', type=float, default=None,
                        help=SUPPRESS)

    ## Store tumor information for debug
    parser.add_argument('--store_tumor_infos', type=str2bool, default=False,
                        help=SUPPRESS)

    ## Output all alternative indel cigar for debug purpose
    parser.add_argument('--alt_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Minimum indel allele frequency for a site to be considered as a candidate site,
    parser.add_argument('--indel_min_af', type=float, default=1.0,
                        help=SUPPRESS)

    ## Minimum INDEL allele frequency for a truth site for training
    parser.add_argument('--min_truth_indel_af', type=float, default=None,
                        help=SUPPRESS)

    ## Truth vcf fn, use for training
    parser.add_argument('--truth_vcf_fn', type=str, default=None,
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

    parser.add_argument('--flanking', type=int, default=None,
                        help=SUPPRESS)


    args = parser.parse_args()

    extract_pair_candidates(args)


if __name__ == "__main__":
    main()
