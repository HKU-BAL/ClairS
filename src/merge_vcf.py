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

import subprocess
import os
from math import log, e
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader, VcfWriter
import shared.param as param
from shared.utils import log_warning, str2bool, str_none

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]
Phred_Trans = (-10 * log(e, 10))

def quality_score_from(probability, int_format=False, use_phred_qual=True):
    p = float(probability)
    if use_phred_qual:
        tmp = max(Phred_Trans * log(((1.0 - p) + 1e-10) / (p + 1e-10)) + 2.0, 0.0)
    else:
        tmp = max(p, 0.0)

    return str(int(round(tmp, 3))) if int_format else "%.3f" % float(round(tmp, 3))

def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

def mark_low_qual(row, quality_score_for_pass):
    if row == '' or "Germline" in row or "RefCall" in row:
        return row
    columns = row.split('\t')
    qual = float(columns[5])
    if quality_score_for_pass and qual <= quality_score_for_pass:
        columns[6] = "LowQual"

    return '\t'.join(columns)


def update_GQ(columns):
    INFO = columns[8]
    FORMAT = columns[9].split(':')
    gq_index = INFO.split(':').index("GQ")
    FORMAT[gq_index] = str(int(float(columns[5]))) if float(columns[5]) > 0.0 else quality_score_from(FORMAT[gq_index], True)
    columns[9] = ':'.join(FORMAT)
    return columns

def merge_vcf(args):
    compress_vcf = args.compress_vcf
    platform = args.platform
    use_phred_qual = args.use_phred_qual
    cmdline_file = args.cmdline
    prefer_recall = args.prefer_recall

    cmdline = None
    if cmdline_file is not None and os.path.exists(cmdline_file):
        cmdline = open(cmdline_file).read().rstrip()

    max_qual_filter_fa_calls = args.max_qual_filter_fa_calls if args.max_qual_filter_fa_calls is not None else param.qual_dict[platform]
    quality_score_for_pass = args.qual if args.qual is not None else param.qual_dict[platform]
    af_cut_off = args.af if args.af is not None else param.af_dict[platform]
    input_vcf_reader = VcfReader(vcf_fn=args.full_alignment_vcf_fn,
                                 ctg_name=None,
                                 show_ref=True,
                                 keep_row_str=True,
                                 skip_genotype=True,
                                 filter_tag=None,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    fa_input_variant_dict = input_vcf_reader.variant_dict

    pass_fa_set = set([k for k, v in fa_input_variant_dict.items() if v.filter == "PASS"])

    pileup_input_variant_dict = defaultdict(str)
    row_count = 0
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    filter_count = 0
    af_filter_count = 0
    f = open(args.pileup_vcf_fn)
    recall_count = 0
    for row in f:
        if row[0] == '#':
            continue

        row_count += 1
        columns = row.strip().split()
        ctg_name, pos = columns[0], columns[1]
        qual = float(columns[5])
        filter = columns[6]
        if filter != 'PASS':
            k = (ctg_name, int(pos))
            pileup_input_variant_dict[k] = row
            continue

        if max_qual_filter_fa_calls is not None:
            if qual <= max_qual_filter_fa_calls:
                if (ctg_name, int(pos)) not in pass_fa_set:
                    if prefer_recall:
                        columns[5] = quality_score_from(columns[5], use_phred_qual=use_phred_qual)
                        columns = update_GQ(columns)
                        row = '\t'.join(columns) + '\n'
                        contig_dict[ctg_name][int(pos)] = row
                        recall_count += 1
                    # added pileup and fa records into output
                    else:
                        filter_count += 1
                    continue
            elif (ctg_name, int(pos)) not in pass_fa_set and platform == 'ilmn':
                columns[5] = quality_score_from(columns[5], use_phred_qual=use_phred_qual)
                columns = update_GQ(columns)
                row = '\t'.join(columns) + '\n'
                contig_dict[ctg_name][int(pos)] = row

        if af_cut_off is not None:
            tag_list = columns[8].split(':')
            taf_index = tag_list.index('AF') if 'AF' in tag_list else tag_list.index('VAF')
            af = float(columns[9].split(':')[taf_index])
            if af <= af_cut_off and (ctg_name, int(pos)) not in pass_fa_set:
                af_filter_count += 1
                continue

        if (ctg_name, int(pos)) in pass_fa_set:
            QUAL = (qual + float(fa_input_variant_dict[(ctg_name, int(pos))].qual)) / 2
            columns[5] = quality_score_from(QUAL, use_phred_qual=use_phred_qual)
            #update GQ to phred
            columns = update_GQ(columns)
            row = '\t'.join(columns) + '\n'

            contig_dict[ctg_name][int(pos)] = row
            no_vcf_output = False

    # append all non_pass fa variant if need to print ref calls
    for k, v in fa_input_variant_dict.items():
        if k[0] not in contig_dict or k[1] not in contig_dict[k[0]]:
            row = v.row_str
            columns = row.strip().split()
            if columns[6] != "Germline" and columns[6] != "RefCall":
                if prefer_recall:
                    columns[5] = quality_score_from(columns[5], use_phred_qual=use_phred_qual)
                    columns = update_GQ(columns)
                    row = '\t'.join(columns) + '\n'
                    contig_dict[k[0]][k[1]] = row
                    no_vcf_output = False
                    recall_count += 1
                    continue
                else:
                    columns[5] = "0.000"
            else:
                columns[5] = quality_score_from(columns[5], use_phred_qual=use_phred_qual)
            #update GQ to phred
            columns = update_GQ(columns)
            row = '\t'.join(columns) + '\n'
            contig_dict[k[0]][k[1]] = row
            no_vcf_output = False

    # append all non_pass pileup variant if need to print ref calls
    for k, row in pileup_input_variant_dict.items():
        if k[0] not in contig_dict or k[1] not in contig_dict[k[0]]:
            columns = row.strip().split()
            if columns[6] != "Germline" and columns[6] != "RefCall":
                columns[5] = "0.000"
            else:
                columns[5] = quality_score_from(columns[5], use_phred_qual=use_phred_qual)
            #update GQ to phred
            columns = update_GQ(columns)
            row = '\t'.join(columns) + '\n'
            contig_dict[k[0]][k[1]] = row
            no_vcf_output = False

    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, please check the setting"))
    if no_vcf_output:
        print(log_warning("[WARNING] No variant found, please check the setting"))

    print("[INFO] Full-alignment variants filtered by pileup: ", filter_count)
    if args.af is not None:
        print("[INFO] Full-alignment variants filtered by AF: ", af_filter_count)
    if prefer_recall:
        print("[INFO] --prefer_recall enabled! Total recalled records: ", recall_count)
    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))

    output_vcf_writer = VcfWriter(vcf_fn=args.output_fn,
                                 ctg_name=','.join(list(contig_dict.keys())),
                                 ref_fn=args.ref_fn,
                                 sample_name=args.sample_name,
                                 cmdline=cmdline,
                                 show_ref_calls=True)

    for contig in contigs_order_list:
        all_pos = sorted(contig_dict[contig].keys())
        for pos in all_pos:
            #Mark low QUAL
            row = mark_low_qual(contig_dict[contig][pos], quality_score_for_pass)
            output_vcf_writer.vcf_writer.write(row)
    output_vcf_writer.close()

    if compress_vcf:
        compress_index_vcf(args.output_fn)

    if args.enable_indel_calling and not args.indel_calling:
        output_fn = args.output_fn + '.gz' if compress_vcf else args.output_fn
        output_dir = os.path.dirname(output_fn)
        file_name = output_fn.split('/')[-1]
        if compress_vcf:
            subprocess.run("cd {} && ln -sf {} snv.vcf.gz".format(output_dir, file_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subprocess.run("cd {} && ln -sf {}.tbi snv.vcf.gz.tbi".format(output_dir, file_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            subprocess.run("cd {} && ln -sf {} snv.vcf.vcf".format(output_dir, file_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def main():
    parser = ArgumentParser(description="Merge full-alignment with pileup VCF")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Output VCF filename, required")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Full-alignment VCF input")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--compress_vcf', type=str2bool, default=True,
                        help="Compress and index the output VCF")

    parser.add_argument('--cmdline', type=str_none, default=None,
                        help="If defined, added command line into VCF header")

    # options for advanced users
    parser.add_argument('--qual', type=float, default=None,
                        help="EXPERIMENTAL: If set, variants Phread quality with >=$qual will be marked 'PASS', or 'LowQual' otherwise")

    parser.add_argument('--use_phred_qual', type=str2bool, default=True,
                        help="EXPERIMENTAL: Use Phred quality score instead of probability")

    parser.add_argument('--enable_indel_calling', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable Indel calling, make the snv output vcf file soft link")

    parser.add_argument('--af', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--max_qual_filter_fa_calls', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--bed_format', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--indel_calling', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--prefer_recall', type=str2bool, default=False,
                        help=SUPPRESS)

    args = parser.parse_args()

    merge_vcf(args)


if __name__ == "__main__":
    main()
