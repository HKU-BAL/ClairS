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

from sys import stderr
from subprocess import Popen
from subprocess import PIPE
from argparse import ArgumentParser
from collections import defaultdict

from shared.interval_tree import bed_tree_from, is_region_in

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
        if bed_tree and not is_region_in(tree=bed_tree, contig_name=contig_name,region_start=pos):
            continue

        genotype_info = columns[9].split(':')
        genotype = genotype_info[0]
        g1, g2 = genotype.replace('|', "/").split('/')
        if g1 == "1" and g2 == "1":
            homo_variant_set.add(pos)
            homo_variant_info[pos] = (ref_base,alt_base)
        variant_set.add(pos)
    return homo_variant_set, homo_variant_info, variant_set

split_bed_size = 2000


def decode_af(input_dir, file_list, output_depth=False, output_alt=False, bed_tree=None, contig_name=None):
    alt_info_dict = defaultdict()
    pos_set = set()
    for f in file_list:
        f = os.path.join(input_dir, f)
        if not os.path.exists(f):
            print('{} not exist'.format(f))
        for row in open(f):
            if output_alt:
                row = row.rstrip().split('\t')
                if len(row) < 5:
                    continue
                pos = row[1]

                alt_infos = dict([item.split(':') for item in row[5].split(' ')])
                alt_info_dict[int(pos)] = alt_infos
            else:
                row = row.rstrip().split()
                pos = row[1]
                depth = row[3]
                if bed_tree and not is_region_in(bed_tree, contig_name, int(pos)):
                    continue
                min_depth = 4
                if int(depth) < min_depth:
                    continue
                pos_set.add(int(pos))
    if output_alt:
        return alt_info_dict
    return pos_set

def find_candidate_match(alt_info_dict, ref_base, alt_base):
    alt_base = alt_base[0]
    #snp
    if len(ref_base) == 1 and len(alt_base) == 1:
        if alt_base in alt_info_dict:
            return alt_info_dict[alt_base], 'snp'

    # insertion
    if len(ref_base) == 1 and len(alt_base) > 1:
        ab = alt_base[0] + '+' + alt_base[1:]
        if ab in alt_info_dict:
            return alt_info_dict[ab], 'ins'
    # deletion
    if len(ref_base) > 1 and len(alt_base) == 1:
        ab = ref_base[0] + '-' + 'N' * len(ref_base[1:])
        if ab in alt_info_dict:
            return alt_info_dict[ab], 'del'

    if len(alt_info_dict) > 0:
        alt_list = sorted(list(alt_info_dict.items()), key=lambda x: x[1], reverse=True)
        return alt_list[0][1], 'signal'
    return None, ""


def find_tumor_truth_in_normal(args):

    contig_name =args.ctg_name
    #tumor bed fn
    bed_fn = args.bed_fn
    normal_sample = args.normal_sample
    tumor_sample = args.tumor_sample
    reference_cans = args.reference_cans
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

    file_list = os.listdir(normal_alt_dir)
    normal_file_list = [f for f in os.listdir(normal_alt_dir) if "_" + contig_name + "_" in f and f.startswith(normal_sample)]
    tumor_file_list = [f for f in os.listdir(tumor_alt_dir) if "_" + contig_name + "_" in f and f.startswith(tumor_sample)]

    bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=contig_name)
    normal_alt_info_dict = decode_af(normal_alt_dir, normal_file_list, output_alt=True, bed_tree=bed_tree, contig_name=contig_name)
    tumor_alt_info_dict = decode_af(tumor_alt_dir, tumor_file_list, output_alt=True, bed_tree=bed_tree, contig_name=contig_name)

    count_dict = {'total':0, 'snp':0, 'ins':0, 'del':0,'signal':0}
    for pos in unified_variant_dict:
        if not is_region_in(bed_tree,contig_name,pos) or pos in normal_unified_variant_dict:
            continue
        count_dict['total'] += 1
        if unified_variant_dict[pos].genotype == [1, 1] and pos in normal_alt_info_dict:
            ref_base, alt_base = unified_variant_dict[pos].reference_bases, unified_variant_dict[pos].alternate_bases
            alt_info_dict = normal_alt_info_dict[pos]
            af, variant_type = find_candidate_match(alt_info_dict, ref_base, alt_base)
            if af is not None:
                count_dict[variant_type] += 1
    print (count_dict.items())

def main():
    parser = ArgumentParser(description="Find tumor truth variants INFO in normal")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input, default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--normal_alt_dir', type=str, default="PIPE",
                        help="Normal alternative directory")

    parser.add_argument('--tumor_alt_dir', type=str, default="PIPE",
                        help="Tumor alternative directory")

    parser.add_argument('--normal_sample', type=str, default=None,
                        help="Normal sample name")

    parser.add_argument('--tumor_sample', type=str, default=None,
                        help="Tumor sample name")

    parser.add_argument('--unified_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input")

    parser.add_argument('--normal_unified_vcf_fn', type=str, default=None,
                        help="Normal candidate sites VCF file input")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--reference_cans', type=str, default=None,
                        help="The path of all reference files")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctg_name and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctg_name and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--add_truths', action='store_true',
                        help="Include all truths variants")


    args = parser.parse_args()

    find_tumor_truth_in_normal(args)


if __name__ == "__main__":
    main()
