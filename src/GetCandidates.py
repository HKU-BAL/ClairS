import shlex
import os
import sys
from os.path import isfile, abspath
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
import argparse
import shlex
from subprocess import PIPE
import subprocess
from collections import defaultdict

from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict
from shared.interval_tree import bed_tree_from, is_region_in
import shared.param as param

from shared.utils import AltInfos
def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def vcf_reader(vcf_fn, contig_name, bed_tree=None, add_hete_pos=False):
    homo_variant_set = set()
    variant_set = set()
    variant_info = defaultdict()
    homo_variant_info = defaultdict()
    hete_variant_set = set()
    hete_variant_info = defaultdict()

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
        if bed_tree and not is_region_in(tree=bed_tree, contig_name=contig_name,region_start=pos):
            continue

        genotype_info = columns[9].split(':')
        genotype = genotype_info[0]
        g1, g2 = genotype.replace('|', "/").split('/')
        if g1 == "1" and g2 == "1":
            homo_variant_set.add(pos)
            homo_variant_info[pos] = (ref_base,alt_base)
        if add_hete_pos and (g1 == "1" and g2 == "0" or g1 == "0" and g2 == "1" or g1 == "1" and g2 == "2"):
            hete_variant_set.add(pos)
            hete_variant_info[pos] = (ref_base, alt_base)
        variant_set.add(pos)
        variant_info[pos] = (ref_base, alt_base)
    return homo_variant_set, homo_variant_info, hete_variant_set, hete_variant_info, variant_set, variant_info

split_bed_size = 10000

def get_ref_candidates(fn, contig_name = None, bed_tree=None):
    ref_cans_dict = defaultdict(AltInfos)
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

        af_list = af_infos.split(',')
        alt_dict = dict([[item.split(':')[0], float(item.split(':')[1])] for item in alt_infos.split(' ')])
        ref_cans_dict[pos] = AltInfos(pos=pos,
                                      ref_base=ref_base,
                                      depth=depth,
                                      af_list=af_list,
                                      alt_dict=alt_dict)
    return ref_cans_dict


def find_candidate_match(alt_info_dict, ref_base, alt_base):
    # alt_base = alt_base[0]
    alt_list = sorted(list(alt_info_dict.items()), key=lambda x: x[1], reverse=True)
    max_af = alt_list[0][1] if len(alt_info_dict) > 0 else None
    #snp
    if len(ref_base) == 1 and len(alt_base) == 1:
        if alt_base in alt_info_dict:
            return alt_info_dict[alt_base], 'snp', max_af
    #
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


def filter_somatic_candidates(truths, variant_info, alt_dict, paired_alt_dict):
    filtered_truths = []
    count_dict = {'total': 0, 'snp': 0, 'ins': 0, 'del': 0, 'signal': 0}
    truth_not_pass_af = 0
    truth_filter_in_normal = 0
    for pos, variant_type in truths:
        if pos not in alt_dict:
            truth_not_pass_af += 1

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
            #     'truth_not_pass_af')
            # print(vcf_format)

            continue

        if pos in paired_alt_dict:
            ref_base, alt_base = variant_info[pos]
            af, vt, max_af = find_candidate_match(alt_info_dict=paired_alt_dict[pos].alt_dict, ref_base=ref_base, alt_base=alt_base)
            if af is not None:
                count_dict[vt] += 1
                truth_filter_in_normal += 1
                continue
            # if max_af > 0.2:
            #     vcf_format = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF:VT\t%s:%d:%d:%.4f:%s" % (
            #         'chr1',
            #         int(pos),
            #         "A",
            #         "A",
            #         10,
            #         'PASS',
            #         '.',
            #         "0/0",
            #         10,
            #         10,
            #         0.5,
            #         'truth_not_pass_af')
            #     print(vcf_format)

        filtered_truths.append([pos, variant_type])
    # print (count_dict, truth_not_pass_af)
    print (truth_filter_in_normal)
    return filtered_truths


def get_candidates(args):
    contig_name = args.ctgName
    bed_fn = args.bed_fn
    bed_tree = bed_tree_from(bed_fn, contig_name)
    vcf_fn_1 = args.vcf_fn_1
    vcf_fn_2 = args.vcf_fn_2
    normal_reference_cans_fn = args.normal_reference_cans
    tumor_reference_cans_fn = args.tumor_reference_cans
    fp_fn = args.fp_fn
    tp_fn = args.tp_fn
    add_hete_pos = args.add_hete_pos
    flankingBaseNum = args.flankingBaseNum if args.flankingBaseNum else param.flankingBaseNum
    split_folder = args.split_folder
    homo_variant_set_1, homo_variant_info_1, hete_variant_set_1, hete_variant_info_1, variant_set_1, variant_info_1 = vcf_reader(vcf_fn=vcf_fn_1, contig_name=contig_name, bed_tree=bed_tree, add_hete_pos=add_hete_pos)
    homo_variant_set_2, homo_variant_info_2, hete_variant_set_2, hete_variant_info_2, variant_set_2, variant_info_2 = vcf_reader(vcf_fn=vcf_fn_2, contig_name=contig_name, bed_tree=bed_tree, add_hete_pos=add_hete_pos)
    normal_alt_dict = get_ref_candidates(fn=normal_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree)
    tumor_alt_dict = get_ref_candidates(fn=tumor_reference_cans_fn, contig_name=contig_name, bed_tree=bed_tree)

    normal_ref_cans_list = [pos for pos in normal_alt_dict if pos not in variant_set_2 and pos not in variant_set_1]
    tumor_ref_cans_list = [pos for pos in tumor_alt_dict if pos not in variant_set_2 and pos not in variant_set_1]
    intersection_pos_set = homo_variant_set_1.intersection(homo_variant_set_2)

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

    homo_germline = [(item, 'homo') for item in list(same_alt_pos_set)]
    hete_germline = [(item, 'hete') for item in hete_list_with_same_repre]
    references = [(item, 'ref') for item in normal_ref_cans_list + tumor_ref_cans_list]
    fp_list = sorted(homo_germline + references + hete_germline, key=lambda x: x[0])


    somatic_set = sorted(list(homo_variant_set_2 - homo_variant_set_1))
    hete_somatic_set = sorted(list(hete_variant_set_2 - variant_set_1))
    homo_somatic = [(item, 'homo_somatic') for item in somatic_set if item not in variant_set_1]
    # skip hete variant here
    hete_somatic = [(item, 'hete_somatic') for item in hete_somatic_set] if add_hete_pos else []

    overall_somatic_truths = len(somatic_set)
    overall_somatic_truths_not_in_pair_truth = len(homo_somatic)
    print(len(homo_variant_set_2), overall_somatic_truths, overall_somatic_truths_not_in_pair_truth)
    homo_somatic = filter_somatic_candidates(truths=homo_somatic, variant_info=variant_info_2, alt_dict=tumor_alt_dict, paired_alt_dict=normal_alt_dict)
    hete_somatic = filter_somatic_candidates(truths=hete_somatic, variant_info=variant_info_2, alt_dict=tumor_alt_dict, paired_alt_dict=normal_alt_dict)
    tp_list = sorted(list(homo_somatic + hete_somatic), key=lambda x: x[0])

    # for pos, variant_type in fp_list + tp_list:
    #
    #     vcf_format = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF:VT\t%s:%d:%d:%.4f:%s" % (
    #         contig_name,
    #         int(pos),
    #         "A",
    #         "A",
    #         10,
    #         'PASS',
    #         '.',
    #         "0/0",
    #         10,
    #         10,
    #         0.5,
    #         variant_type)
    #     print(vcf_format)

    for pos_list, file_name in zip([fp_list, tp_list], [fp_fn, tp_fn]):
        all_full_aln_regions = []
        region_num = len(pos_list) // split_bed_size + 1 if len(
            pos_list) % split_bed_size else len(pos_list) // split_bed_size

        for idx in range(region_num):
            # a windows region for create tensor # samtools mpileup not include last position
            split_output = pos_list[idx * split_bed_size: (idx + 1) * split_bed_size]

            split_output = [(item[0] - flankingBaseNum, item[0] + flankingBaseNum + 2, item[1]) for item in
                            split_output]

            output_path = os.path.join(split_folder, file_name, '{}.{}_{}_{}'.format(contig_name, idx, region_num, file_name))
            all_full_aln_regions.append(output_path)
            with open(output_path, 'w') as output_file:
                output_file.write('\n'.join(
                    ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), x[2]]) for x in
                     split_output]) + '\n')  # bed format

        all_full_aln_regions_path = os.path.join(split_folder, file_name, 'FULL_ALN_FILE_{}'.format( contig_name))
        with open(all_full_aln_regions_path, 'w') as output_file:
            output_file.write('\n'.join(all_full_aln_regions) + '\n')

    print_all = False
    if print_all:
        print ('[INFO] {} total homo reference pos: 1:{}, 2:{}, references:{} hete variants:{} homo truth with same pos intersection:{}, homo truth with_same_repre:{}'.format(contig_name, len(homo_variant_set_1), len(homo_variant_set_2), len(ref_cans_list), len(hete_list_with_same_repre), len(intersection_pos_set), len(same_alt_pos_set)))
    else:
        print ('[INFO] {} homo germline:{} hete germline:{} references:{} homo somatic:{} hete somatic:{}'.format(contig_name, len(homo_germline), len(hete_germline), len(references),len(homo_somatic),len(hete_somatic) ))


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",  # required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",  # required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn_1', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--vcf_fn_2', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.08,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--normal_reference_cans', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--tumor_reference_cans', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--fp_fn', type=str, default="fp",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--tp_fn', type=str, default="tp",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

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

    parser.add_argument('--flankingBaseNum', type=int, default=None,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    # options for debug purpose
    parser.add_argument('--add_hete_pos', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")


    args = parser.parse_args()

    get_candidates(args)


if __name__ == "__main__":
    main()
