import os
import subprocess
import shlex
from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict
import sys
import shlex

from shared.vcf import VcfReader
import shared.param as param
from shared.utils import log_error, log_warning, file_path_from, str2bool

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def merge_vcf(args):
    compress_vcf = args.compress_vcf
    platform = args.platform
    qual_cut_off = args.qual if args.qual is not None else param.qual_dict[platform]
    af_cut_off = args.af if args.af is not None else param.af_dict[platform]

    input_vcf_reader = VcfReader(vcf_fn=args.full_alignment_vcf_fn, ctg_name=None, show_ref=False, keep_row_str=True,
                                 skip_genotype=True,
                                 filter_tag=None,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    fa_input_variant_dict = input_vcf_reader.variant_dict

    pass_fa_set = set([k for k, v in fa_input_variant_dict.items() if v.filter == "PASS"])

    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    filter_count = 0
    af_filter_count = 0
    f = open(args.pileup_vcf_fn)
    for row in f:

        if row[0] == '#':
            if row not in header:
                header.append(row)
            continue
        row_count += 1
        # use the first vcf header
        columns = row.strip().split()
        ctg_name, pos = columns[0], columns[1]
        qual = float(columns[5])
        filter = columns[6]
        if filter != 'PASS':
            continue

        if qual_cut_off is not None and qual <= qual_cut_off:
            if (ctg_name, int(pos)) not in pass_fa_set:
                filter_count += 1
                continue

        if af_cut_off is not None:
            tag_list = columns[8].split(':')
            taf_index = tag_list.index('AF') if 'AF' in tag_list else tag_list.index('VAF')
            af = float(columns[9].split(':')[taf_index])
            if af <= af_cut_off and (ctg_name, int(pos)) not in pass_fa_set:
                af_filter_count += 1
                continue

        if (ctg_name, int(pos)) in pass_fa_set:
            columns[5] = str((qual + float(fa_input_variant_dict[(ctg_name, int(pos))].qual)) / 2)
            row = '\t'.join(columns) + '\n'

        contig_dict[ctg_name][int(pos)] = row
        no_vcf_output = False

    # append all non_pass variant if need to print ref calls
    for k, v in fa_input_variant_dict.items():
        if k[0] not in contig_dict or k[1] not in contig_dict[k[0]]:
            contig_dict[k[0]][k[1]] = v.row_str

    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, please check the setting"))
    if no_vcf_output:
        print(log_warning("[WARNING] No variant found, please check the setting"))

    print("Filter count", filter_count)
    if args.af is not None:
        print("AF filter count", af_filter_count)
    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(args.output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])
    f.close()

    if compress_vcf:
        compress_index_vcf(args.output_fn)


def main():
    parser = ArgumentParser(description="Sort a VCF file according to contig name and starting position")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input, default: %(default)s")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Output VCF filename, required")

    parser.add_argument('--input_dir', type=str, default=None,
                        help="Input directory")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Full-alignment VCF input")

    parser.add_argument('--vcf_fn_prefix', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--vcf_fn_suffix', type=str, default='.vcf',
                        help="Input vcf filename suffix")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--contigs_fn', type=str, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--bed_format', action='store_true',
                        help="Only work for gvcf file, reduce hard disk space")

    parser.add_argument('--filter_vcf_fn', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--qual', type=float, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--af', type=float, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--compress_vcf', type=str2bool, default=True,
                        help="Compress and index output VCF")

    args = parser.parse_args()

    merge_vcf(args)


if __name__ == "__main__":
    main()
