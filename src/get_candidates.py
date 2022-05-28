import shlex
import os
import sys
import math
from os.path import isfile, abspath
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
import argparse
import random
import shlex
import math
from subprocess import PIPE
from shared.vcf import VcfWriter
from shared.bed import BedWriter
import subprocess
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
        if bed_tree and not is_region_in(tree=bed_tree, contig_name=contig_name,region_start=pos):
            continue

        genotype_info = columns[9].split(':')
        genotype = genotype_info[0]
        g1, g2 = genotype.replace('|', "/").split('/')
        if g1 == "1" and g2 == "1":
            homo_variant_set.add(pos)
            homo_variant_info[pos] = (ref_base, alt_base)
        if add_hetero_pos and (g1 == "1" and g2 == "0" or g1 == "0" and g2 == "1"): #or g1 == "1" and g2 == "2"): # exclude 1/2 truths here
            hetero_variant_set.add(pos)
            hetero_variant_info[pos] = (ref_base, alt_base)

    return homo_variant_set, homo_variant_info, hetero_variant_set, hetero_variant_info, variant_set, variant_info

def get_ref_candidates(fn, contig_name = None, bed_tree=None, variant_info=None):
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

            # if variant_info is not None and pos not in variant_info:
            #     continue
            ref_base, depth, af_infos, alt_infos = columns[2:6]
            tumor_infos = columns[6] if len(columns) > 6 else ""

            af_list = af_infos.split(',')
            alt_dict = dict([[item.split(':')[0], float(item.split(':')[1])] for item in alt_infos.split(' ')])
            tumor_alt_dict = dict([[item.split(':')[0], float(item.split(':')[1])] for item in tumor_infos.split(' ')]) if len(tumor_infos) else dict()
            ref_cans_dict[pos] = AltInfos(pos=pos,
                                          ref_base=ref_base,
                                          depth=depth,
                                          af_list=af_list,
                                          alt_dict=alt_dict,
                                          tumor_alt_dict=tumor_alt_dict)
        unzip_process.stdout.close()
        unzip_process.wait()
    return ref_cans_dict


def find_most_frequent_candidate(alt_info_dict, ref_base):
    alt_list = sorted(list([k,v] for k,v in alt_info_dict.items() if k != ref_base), key=lambda x: x[1], reverse=True)

    for alt_base, af in alt_list:
        if alt_base == ref_base or len(alt_base) > 1:
            continue
        return alt_base, af

    return None, None

def find_candidate_match(alt_info_dict, ref_base, alt_base):
    # alt_base = alt_base[0]
    alt_list = sorted(list(alt_info_dict.items()), key=lambda x: x[1], reverse=True)
    max_af = alt_list[0][1] if len(alt_info_dict) > 0 else None
    #snp
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


def filter_germline_candidates(truths, variant_info, alt_dict, paired_alt_dict, gen_vcf, INFO):
    #alt_dict: tumor
    #pair alt_dict: normal
    filtered_truths = []
    low_confident_germline_truths = []
    truth_not_pass_af = 0
    germline_filtered_by_af_distance = 0
    for pos, variant_type in truths:
        if pos not in alt_dict:
            truth_not_pass_af += 1
            if gen_vcf:
                low_confident_germline_truths.append((pos, variant_type + INFO + 'germline_no_read_support'))
            # vcf_format = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF:VT\t%s:%d:%d:%.4f:%s" % (
            #     'chr1',
            #     int(pos),
            #     "A",
            #     "A",
            #     10,
            #     'PASS',
            #     '.',
            #     "0/0",
            #     10,
            #     10,
            #     0.5,
            #     variant_type)
            # print(vcf_format)
            continue

        # discard the germline variants, if normal and tumor have large distinct
        if pos in paired_alt_dict:
            ref_base, alt_base = variant_info[pos]
            pair_af, vt, max_af = find_candidate_match(alt_info_dict=paired_alt_dict[pos].alt_dict, ref_base=ref_base, alt_base=alt_base)
            paired_alt_dict[pos].max_candidate_af = max_af
            paired_alt_dict[pos].support_alternative_af = pair_af

            af, vt, max_af = find_candidate_match(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base, alt_base=alt_base)
            alt_dict[pos].max_candidate_af = max_af
            alt_dict[pos].support_alternative_af = pair_af
            # pair_af: af in normal, af: af in tumor
            if pair_af is None or af is None or math.fabs(pair_af - af) > 0.1:
                germline_filtered_by_af_distance += 1
                if gen_vcf:
                    low_confident_germline_truths.append((pos, variant_type + INFO + 'germline_af_gap'))
                continue

            if pos in paired_alt_dict and paired_alt_dict[pos].support_alternative_af is None:
                paired_alt_dict[pos].support_alternative_af = pair_af
            if pos in alt_dict and alt_dict[pos].support_alternative_af is None:
                alt_dict[pos].support_alternative_af = af
            filtered_truths.append([pos, variant_type])
    # print (truth_not_pass_af, germline_filtered_by_af_distance)

    print(
        '[INFO] {} truth germline variants filtered by no read support in normal BAM: {}, filtered by large AF distance between normal and tumor: {}'.format(
            INFO, truth_not_pass_af, germline_filtered_by_af_distance))


    return filtered_truths, low_confident_germline_truths



def filter_reference_candidates(truths, alt_dict, paired_alt_dict, gen_vcf, INFO):
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
        most_frequent_alt_base, tumor_af = find_most_frequent_candidate(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base)


        if most_frequent_alt_base is None or tumor_af is None:
            continue
        tumor_af = float(tumor_af)
        if pos not in paired_alt_dict:
            filtered_truths.append([pos, variant_type])
        else:
            pair_af, vt, pair_max_af = find_candidate_match(alt_info_dict=paired_alt_dict[pos].alt_dict, ref_base=ref_base, alt_base=most_frequent_alt_base)
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



def filter_somatic_candidates(truths, variant_info, alt_dict, paired_alt_dict, gen_vcf=False, INFO="Homo"):
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
        if pos not in alt_dict:# or not len(alt_dict[pos].tumor_alt_dict):
            truth_not_pass_af += 1
            if gen_vcf:
                low_confident_truths.append((pos, variant_type + INFO + '_no_read_support'))
            if skip_no_read_support:
                continue

        normal_af, tumor_reads_af = None, None

        if pos in paired_alt_dict:
            ref_base, alt_base = variant_info[pos]
            # very high af with same alt_base in normal is not suitable as candidate for training
            normal_af, vt, max_af = find_candidate_match(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base, alt_base=alt_base)
            if (normal_af is not None and normal_af > max_af_in_normal) or (max_af is not None and max_af > max_af_in_normal):
                truth_filter_in_normal += 1
                if gen_vcf:
                    low_confident_truths.append((pos, variant_type + INFO + '_high_af_in_normal'))
                if skip_high_af_in_normal:
                    continue


        if len(alt_dict[pos].tumor_alt_dict):
            ref_base, alt_base = variant_info[pos]
            # tumor_af, _, _ = find_candidate_match(alt_info_dict=alt_dict[pos].alt_dict, ref_base=ref_base,alt_base=alt_base)
            tumor_reads_af, vt, max_af = find_candidate_match(alt_info_dict=alt_dict[pos].tumor_alt_dict, ref_base=ref_base, alt_base=alt_base)
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
        # if math.fabs(tumor_af - tumor_af_from_tumor_reads - normal_af) > af_gap_for_errors:
        #     continue

        if pos in paired_alt_dict and paired_alt_dict[pos].support_alternative_af is None:
            paired_alt_dict[pos].support_alternative_af = normal_af
        if pos in alt_dict and alt_dict[pos].support_alternative_af is None:
            alt_dict[pos].support_alternative_af = tumor_reads_af

        filtered_truths.append([pos, variant_type])

    print ('[INFO] {} truth variants filtered by high AF in normal BAM: {}, filtered by low AF {}:{}, filtered by no/low read support:{}/{}'.format(INFO, truth_filter_in_normal, min_af_for_tumor, truth_filter_with_low_af, truth_not_pass_af, truth_filter_with_low_coverage))
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
    proportion = args.proportion
    synthetic_coverage = args.synthetic_coverage
    output_bed_fn = args.output_bed_fn
    output_fp_bed_regions = output_bed_fn is not None
    homo_variant_set_1, homo_variant_info_1, hete_variant_set_1, hete_variant_info_1, variant_set_1, variant_info_1 = vcf_reader(vcf_fn=vcf_fn_1, contig_name=contig_name, bed_tree=bed_tree, add_hete_pos=add_hete_pos)
    homo_variant_set_2, homo_variant_info_2, hete_variant_set_2, hete_variant_info_2, variant_set_2, variant_info_2 = vcf_reader(vcf_fn=vcf_fn_2, contig_name=contig_name, bed_tree=bed_tree, add_hete_pos=add_hete_pos)
    # print (len(variant_info_1), len(variant_info_2))
    tumor_alt_dict = get_ref_candidates(fn=tumor_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree, variant_info=variant_info_2)
    normal_alt_dict = get_ref_candidates(fn=normal_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree, variant_info=variant_info_2)

    normal_ref_cans_list = [pos for pos in normal_alt_dict if pos not in variant_set_2 and pos not in variant_set_1]
    tumor_ref_cans_list = [pos for pos in tumor_alt_dict if pos not in variant_set_2 and pos not in variant_set_1]
    intersection_pos_set = homo_variant_set_1.intersection(homo_variant_set_2)

    #
    same_alt_pos_set = set()
    for pos in intersection_pos_set:
        if homo_variant_info_1[pos] != homo_variant_info_2[pos]:
            continue
        same_alt_pos_set.add(pos)

    hete_list_with_same_repre = []
    if add_hete_pos:
        for pos, (ref_base, alt_base) in hete_variant_info_2.items():
            if pos in hete_variant_set_1:
                ref_base_2, alt_base_2 = hete_variant_info_1[pos]
                if ref_base == ref_base_2 and alt_base == alt_base_2:
                    hete_list_with_same_repre.append(pos)

    homo_germline = [(item, 'homo_germline') for item in list(same_alt_pos_set)]
    hete_germline = [(item, 'hete_germline') for item in hete_list_with_same_repre]
    ref_list = normal_ref_cans_list + tumor_ref_cans_list if consider_normal_af else tumor_ref_cans_list
    references = [(item, 'ref') for item in ref_list]
    homo_germline = filter_germline_candidates(truths=homo_germline, variant_info=variant_info_2, alt_dict=tumor_alt_dict, paired_alt_dict=normal_alt_dict, INFO="Homo")
    # need add hete, otherwise, in real case, the performance of hetero variants are too bad
    hete_germline = filter_germline_candidates(truths=hete_germline, variant_info=variant_info_2, alt_dict=tumor_alt_dict, paired_alt_dict=normal_alt_dict, INFO="Hetero")

    add_germline = True
    if not add_germline or gen_vcf:
        # exclude germline variants if we exclude germline variants into training
        # exclude gerlmine variants when create VCF and fp BED region
        homo_germline = []
        hete_germline = []
    if maximum_non_variant_ratio is not None:
        # random sample reference calls in training mode, with seed
        random.seed(0)
        references = random.sample(references, int(len(references) * maximum_non_variant_ratio))
    fp_list = homo_germline + references + hete_germline

    homo_somatic_set = sorted(list(homo_variant_set_2 - variant_set_1))
    hete_somatic_set = sorted(list(hete_variant_set_2 - variant_set_1))
    homo_somatic = [(item, 'homo_somatic') for item in homo_somatic_set]
    # skip hete variant here
    hete_somatic = [(item, 'hete_somatic') for item in hete_somatic_set] if add_hete_pos else []

    homo_somatic, homo_low_confident_truths = filter_somatic_candidates(truths=homo_somatic,
                                             variant_info=variant_info_2,
                                             alt_dict=tumor_alt_dict,
                                             paired_alt_dict=normal_alt_dict,
                                             gen_vcf=gen_vcf)
    hete_somatic, hete_low_confident_truths = filter_somatic_candidates(truths=hete_somatic,
                                             variant_info=variant_info_2,
                                             alt_dict=tumor_alt_dict,
                                             paired_alt_dict=normal_alt_dict,
                                             gen_vcf=gen_vcf,
                                             INFO="Hete")

    # hete_somatic = filter_somatic_candidates(truths=hete_somatic, variant_info=variant_info_2, alt_dict=tumor_alt_dict, paired_alt_dict=normal_alt_dict)
    tp_list = homo_somatic + hete_somatic
    pos_list = sorted(fp_list + tp_list, key=lambda x: x[0])

    if gen_vcf:
        vcf_writer = VcfWriter(vcf_fn=output_vcf_fn, ref_fn=ref_fn, ctg_name=contig_name, show_ref_calls=True)
        for pos, variant_type in pos_list + homo_low_confident_truths + hete_low_confident_truths:
            genotype = '1/1' if (variant_type == 'homo_somatic' or variant_type == 'hete_somatic') else '0/0'
            filter_tag = "PASS" if genotype == '1/1' else "LowQual"
            if genotype == "1/1":
                ref_base, alt_base = variant_info_2[pos]
            elif pos in tumor_alt_dict:
                ref_base = alt_base = tumor_alt_dict[pos].ref_base
            else:
                alt_base = ref_base = normal_alt_dict[pos].ref_base if pos in normal_alt_dict else (tumor_alt_dict[pos] if pos in tumor_alt_dict else "N")

            vcf_writer.write_row(POS=pos,
                                 REF=ref_base,
                                 ALT=alt_base,
                                 QUAL=10,
                                 FILTER=filter_tag,
                                 GT=genotype,
                                 DP=10,
                                 AF=0.5,
                                 VT=variant_type)
        vcf_writer.close()

    if output_fp_bed_regions:
        bed_writer = BedWriter(bed_fn=output_bed_fn)
        for pos, _ in pos_list:
            bed_writer.write_row(ctg_name=contig_name,
                                 start_pos=pos-1,
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

        output_path = os.path.join(split_folder, '{}.{}_{}_por{}_cov{}'.format(contig_name, idx, region_num, int(proportion*100), int(synthetic_coverage)))
        all_candidates_regions.append(output_path + ' ' + str(proportion) + ' ' + str(synthetic_coverage))
        with open(output_path, 'w') as output_file:
            output_file.write('\n'.join(
                ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), x[2]]) for x in
                 split_output]) + '\n')  # bed format

    all_candidate_path = os.path.join(split_folder, 'CANDIDATES_FILE_{}_{}_{}'.format(contig_name, proportion, synthetic_coverage))
    with open(all_candidate_path, 'w') as output_file:
        output_file.write('\n'.join(all_candidates_regions) + '\n')

    # print_all = False
    # if print_all:
    #     print ('[INFO] {} total homo reference pos: 1:{}, 2:{}, references:{} hetero variants:{} homo truth with same pos intersection:{}, homo truth with_same_repre:{}'.format(contig_name, len(normal_homo_variant_set), len(tumor_homo_variant_set), len(ref_cans_list), len(hetero_list_with_same_repre), len(intersection_pos_set), len(same_alt_pos_set)))
    # else:
    print ('[INFO] {}-{} homo germline:{} hetero germline:{} references:{} homo somatic:{} hetero somatic:{}\n'.format(contig_name, proportion, len(homo_germline), len(hetero_germline), len(references),len(homo_somatic),len(hetero_somatic) ))


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",  # required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default=None,  # required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--normal_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--tumor_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.08,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--normal_reference_cans', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--tumor_reference_cans', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctg_start', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctg_name and/or (--ctg_start, --ctgEnd) are set")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctg_name and/or (--ctg_start, --ctgEnd) are set")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # parser.add_argument('--fp_fn', type=str, default="fp",
    #                     help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")
    #
    # parser.add_argument('--tp_fn', type=str, default="tp",
    #                     help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--output_dir', type=str, default="",
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

    parser.add_argument('--flanking_base_num', type=int, default=None,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--sample_normal_af', type=float, default=None,
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--proportion', type=float, default=None,
                        help="DEBUG: The window size for read phasing")

    parser.add_argument('--synthetic_coverage', type=int, default=None,
                        help="DEBUG: The window size for read phasing")


    # options for debug purpose
    parser.add_argument('--add_hetero_pos', type=str2bool, default=0,
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--output_bed_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    ## Maximum non-variant ratio against variant in the training data
    parser.add_argument('--maximum_non_variant_ratio', type=float, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    get_candidates(args)


if __name__ == "__main__":
    main()
