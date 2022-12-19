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


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def vcf_reader(vcf_fn, contig_name, bed_tree=None):
    homo_variant_set = set()
    variant_set = set()
    homo_variant_info = defaultdict()
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
        if bed_tree and not is_region_in(tree=bed_tree, contig_name=contig_name, region_start=pos):
            continue

        genotype_info = columns[9].split(':')
        genotype = genotype_info[0]
        g1, g2 = genotype.replace('|', "/").split('/')
        if g1 == "1" and g2 == "1":
            homo_variant_set.add(pos)
            homo_variant_info[pos] = (ref_base, alt_base)
        variant_set.add(pos)
    return homo_variant_set, homo_variant_info, variant_set


split_bed_size = 2000


def decode_af(input_dir, file_list, output_depth=False, output_alt=False, bed_tree=None, contig_name=None):
    af_dict = defaultdict(float)
    alt_info_dict = defaultdict()
    pos_set = set()
    for f in file_list:
        f = os.path.join(input_dir, f)
        if not os.path.exists(f):
            print('{} not exist'.format(f))
        for row in open(f):
            if output_alt:
                columns = row.rstrip().split('\t')
                if len(columns) < 5:
                    continue
                pos = columns[1]

                af = columns[4].split(',')[0]
                alt_info_dict[int(pos)] = row
                af_dict[int(pos)] = float(af)
            else:
                columns = row.rstrip().split()
                pos = columns[1]
                depth = columns[3]
                if bed_tree and not is_region_in(bed_tree, contig_name, int(pos)):
                    continue
                min_depth = 4
                if int(depth) < min_depth:
                    continue
                pos_set.add(int(pos))
    if output_alt:
        return alt_info_dict
    return pos_set


def filter_ref(args):
    contig_name = args.ctgName
    # tumor bed fn
    bed_fn = args.bed_fn
    normal_sample = args.normal_sample
    tumor_sample = args.tumor_sample
    reference_cans_dir = args.reference_cans_dir
    normal_alt_dir = args.normal_alt_dir
    tumor_alt_dir = args.tumor_alt_dir
    unified_vcf_fn = args.unified_vcf_fn
    add_truths = args.add_truths
    from shared.vcf import VcfReader
    unified_vcf_reader = VcfReader(vcf_fn=unified_vcf_fn, ctg_name=contig_name, is_var_format=False)
    unified_vcf_reader.read_vcf()
    unified_variant_dict = unified_vcf_reader.variant_dict

    normal_unified_vcf_fn = args.normal_unified_vcf_fn
    normal_unified_vcf_reader = VcfReader(vcf_fn=normal_unified_vcf_fn, ctg_name=contig_name, is_var_format=False)
    normal_unified_vcf_reader.read_vcf()
    normal_unified_variant_dict = normal_unified_vcf_reader.variant_dict

    # file_list = os.listdir(normal_alt_dir)
    normal_file_list = [f for f in os.listdir(normal_alt_dir) if
                        "_" + contig_name + "_" in f and f.startswith(normal_sample)]
    tumor_file_list = [f for f in os.listdir(tumor_alt_dir) if
                       "_" + contig_name + "_" in f and f.startswith(tumor_sample)]

    bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=contig_name)
    normal_alt_info_dict = decode_af(normal_alt_dir, normal_file_list, output_alt=True, bed_tree=bed_tree,
                                     contig_name=contig_name)
    tumor_alt_info_dict = decode_af(tumor_alt_dir, tumor_file_list, output_alt=True, bed_tree=bed_tree,
                                    contig_name=contig_name)

    match_count = 0
    normal_reference_cans_fp = open(os.path.join(reference_cans_dir, normal_sample + '_' + contig_name), 'w')
    tumor_reference_cans_fp = open(os.path.join(reference_cans_dir, tumor_sample + '_' + contig_name), 'w')
    pos_in_normal_truth = 0
    pos_in_tumor_truth = 0
    pos_in_normal_truth_and_tumor_truth = 0
    normal_reference_cans_fp.write("#{}\n".format(normal_sample))
    tumor_reference_cans_fp.write("#{}\n".format(tumor_sample))
    for pos, alt_infos in sorted(list(tumor_alt_info_dict.items()), key=lambda x: x[0]):
        if pos in normal_unified_variant_dict:
            pos_in_normal_truth += 1
            if pos in unified_variant_dict:
                pos_in_normal_truth_and_tumor_truth += 1
        # skip position in truth if not add truths
        if pos in unified_variant_dict:
            pos_in_tumor_truth += 1
        if pos in unified_variant_dict and not add_truths:
            continue
        tumor_reference_cans_fp.write(alt_infos)
        match_count += 1
    for pos, alt_infos in sorted(list(normal_alt_info_dict.items()), key=lambda x: x[0]):
        normal_reference_cans_fp.write(alt_infos)
    normal_reference_cans_fp.close()
    tumor_reference_cans_fp.close()

    print(
        '[INFO] {} normal pos/tumor pos/matched: {}/{}/{}, pos in normal truth/pos in tumor truth/pos in normal and truth:{}/{}/{}'.format(
            contig_name, len(normal_alt_info_dict), len(tumor_alt_info_dict), match_count, pos_in_normal_truth,
            pos_in_tumor_truth, pos_in_normal_truth_and_tumor_truth))


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--normal_alt_dir', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--tumor_alt_dir', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--normal_sample', type=str, default=None,
                        help="Normal sample name")

    parser.add_argument('--tumor_sample', type=str, default=None,
                        help="Tumor sample name")

    parser.add_argument('--unified_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--normal_unified_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.08,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--reference_cans_dir', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--output_dir', type=str, default="",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--add_truths', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    args = parser.parse_args()

    filter_ref(args)


if __name__ == "__main__":
    main()
