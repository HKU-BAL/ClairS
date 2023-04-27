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

import os
import subprocess
import shlex
import argparse

from sys import stderr, exit
from subprocess import PIPE
from argparse import ArgumentParser
from collections import defaultdict
from subprocess import Popen

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23))] + [str(a) for a in
                                                                                   list(range(1, 23))]

def str_none(v):
    if v is None:
        return None
    if v.upper() == "NONE":
        return None
    if isinstance(v, str):
        return v

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'ture', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'flase', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


class Position(object):
    def __init__(self,
                 ctg_name=None,
                 pos=None,
                 row_str=None):
        self.ctg_name = ctg_name
        self.pos = pos
        self.row_str = row_str

class AltInfo(object):
    def __init__(self, ref_base='', normal_alt_info="", tumor_alt_info=""):
        self.ref_base = ref_base
        self.normal_alt_info = normal_alt_info
        self.tumor_alt_info = tumor_alt_info

class VcfReader(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 direct_open=False,
                 keep_row_str=False,
                 save_header=False,
                 ):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.variant_dict = defaultdict(Position)
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.header = ""
        self.save_header = save_header

    def read_vcf(self):

        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        if self.direct_open:
            vcf_fp = open(self.vcf_fn)
            vcf_fo = vcf_fp
        else:
            vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
            vcf_fo = vcf_fp.stdout
        for row in vcf_fo:
            columns = row.strip().split()
            if columns[0][0] == "#":
                if self.save_header:
                    self.header += row
                continue

            # position in vcf is 1-based
            chromosome, position = columns[0], int(columns[1])
            key = (chromosome, position) if self.ctg_name is None else position
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue

            row_str = row if self.keep_row_str else False
            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              row_str=row_str)


def get_alt_info(alt_info):

    if alt_info is None or alt_info == "":
        return 0, 0, 0, 0, 0
    try:
        acgt_count = [0, 0, 0, 0]
        depth, seqs = alt_info.split('-')
        seqs = seqs.split(' ')
        alt_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
        for idx, base in enumerate('ACGT'):
            acgt_count[idx] = alt_dict[base] if base in alt_dict else 0
        return int(depth), acgt_count[0], acgt_count[1], acgt_count[2], acgt_count[3]
    except:
        return 0, 0, 0, 0, 0


def switch_genotype_row(row_str, ref_base=None, normal_alt_info=None, tumor_alt_info=None):
    NDP, NAU, NCU, NGU, NTU = get_alt_info(normal_alt_info)
    DP, AU, CU, GU, TU = get_alt_info(tumor_alt_info)
    columns = row_str.rstrip().split('\t')
    if len(columns) < 10:
        columns += ['.'] * (10 - len(columns))
    columns[3] = columns[3][0] if len(columns[3]) > 0 else '.'  # Only keep the reference base for REF
    if ref_base is not None and ref_base != "" and columns[3] != ref_base:
        columns[3] = ref_base
    columns[4] = '.'  # ALT to 0
    columns[5] = "."  # QUAL to .
    columns[6] = "."  # FILTER to .
    columns[7] = "."  # INFO to .
    columns[8] = "GT:DP:NDP:AU:CU:GU:TU:NAU:NCU:NGU:NTU"
    columns[9] = './.' + ":%d:%d" % (DP, NDP) + ":%d:%d:%d:%d" % (AU, CU, GU, TU) + \
                 ":%d:%d:%d:%d" % (NAU, NCU, NGU, NTU)
    row_str = '\t'.join(columns) + '\n'
    return row_str

def decode_hybrid_info(candidates_folder):
    file_list = [item for item in os.listdir(candidates_folder) if item.endswith("hybrid_info")]
    hybrid_info_dict = defaultdict(AltInfo)
    for file in file_list:
        with open(os.path.join(candidates_folder, file)) as f:
            for row in f:
                columns = row.rstrip('\n').split('\t') # ctg_name, pos, ref_base, normal_info, tumor_info
                if len(columns) < 5:
                    continue
                ctg_name, pos, ref_base, normal_alt_info, tumor_alt_info = columns[:5]
                key = (ctg_name, int(pos))
                hybrid_info_dict[key] = AltInfo(ref_base=ref_base,
                                                normal_alt_info=normal_alt_info,
                                                tumor_alt_info=tumor_alt_info)
    return hybrid_info_dict

def genotype_vcf(args):
    genotyping_mode_vcf_fn = args.genotyping_mode_vcf_fn
    hybrid_mode_vcf_fn = args.hybrid_mode_vcf_fn
    call_fn = args.call_fn
    output_fn = args.output_fn
    switch_genotype = args.switch_genotype
    candidates_folder = args.candidates_folder
    if genotyping_mode_vcf_fn is None and hybrid_mode_vcf_fn is None:
        exit("[ERROR] Both VCF {} and additional VCF is None, exit!".format(genotyping_mode_vcf_fn, hybrid_mode_vcf_fn))

    output = open(output_fn, 'w')

    rc = subprocess.run('cp {} {}.bak'.format(call_fn, call_fn), shell=True)

    somatic_vcf_reader = VcfReader(vcf_fn=call_fn,
                                  ctg_name=None,
                                  keep_row_str=True,
                                  save_header=True)

    somatic_vcf_reader.read_vcf()
    somatic_variant_dict = somatic_vcf_reader.variant_dict
    output.write(somatic_vcf_reader.header)

    hybrid_info_dict = defaultdict()
    if candidates_folder is not None and os.path.exists(candidates_folder):
        hybrid_info_dict = decode_hybrid_info(candidates_folder)

    if genotyping_mode_vcf_fn is not None:
        vcf_reader = VcfReader(vcf_fn=genotyping_mode_vcf_fn,
                               ctg_name=None,
                               keep_row_str=True,
                               save_header=True)

        vcf_reader.read_vcf()
        variant_dict = vcf_reader.variant_dict

        all_contigs_list = list(set([k[0] for k in variant_dict]))

        contigs_order = major_contigs_order + all_contigs_list

        contigs_order_list = sorted(all_contigs_list, key=lambda x: contigs_order.index(x))

        count = 0
        contig_dict = defaultdict(list)
        for k, v in variant_dict.items():
            ctg, pos = k
            if k not in somatic_variant_dict:
                row_str = variant_dict[k].row_str
                count += 1
                if switch_genotype:
                    ref_base = None
                    normal_alt_info = None
                    tumor_alt_info = None
                    if k in hybrid_info_dict:
                        ref_base = hybrid_info_dict[k].ref_base
                        normal_alt_info = hybrid_info_dict[k].normal_alt_info
                        tumor_alt_info = hybrid_info_dict[k].tumor_alt_info

                    row_str = switch_genotype_row(row_str,
                                                  ref_base=ref_base,
                                                  normal_alt_info=normal_alt_info,
                                                  tumor_alt_info=tumor_alt_info,)
            else:
                row_str = somatic_variant_dict[k].row_str

            contig_dict[ctg].append((int(pos), row_str))

        for contig in contigs_order_list:
            row_list = [item[1] for item in sorted(contig_dict[contig], key=lambda x: x[0])]
            output.write(''.join(row_list))


    elif hybrid_mode_vcf_fn is not None:
        vcf_reader = VcfReader(vcf_fn=hybrid_mode_vcf_fn,
                               ctg_name=None,
                               keep_row_str=True,
                               save_header=False)

        vcf_reader.read_vcf()
        variant_dict = vcf_reader.variant_dict

        all_contigs_list = list(set([k[0] for k in somatic_variant_dict] + [k[0] for k in variant_dict]))

        contigs_order = major_contigs_order + all_contigs_list

        contigs_order_list = sorted(all_contigs_list, key=lambda x: contigs_order.index(x))

        count = 0
        contig_dict = defaultdict(list)

        #add all called rows first
        for k,v in somatic_variant_dict.items():
            ctg, pos = k
            row_str = somatic_variant_dict[k].row_str
            contig_dict[ctg].append((int(pos), row_str))

        #append the rows not in called records
        for k, v in variant_dict.items():
            ctg, pos = k
            if k not in somatic_variant_dict:
                row_str = variant_dict[k].row_str
                count += 1
                if switch_genotype:
                    ref_base = None
                    normal_alt_info = None
                    tumor_alt_info = None
                    if k in hybrid_info_dict:
                        ref_base = hybrid_info_dict[k].ref_base
                        normal_alt_info = hybrid_info_dict[k].normal_alt_info
                        tumor_alt_info = hybrid_info_dict[k].tumor_alt_info

                    row_str = switch_genotype_row(row_str,
                                                  ref_base=ref_base,
                                                  normal_alt_info=normal_alt_info,
                                                  tumor_alt_info=tumor_alt_info)

                contig_dict[ctg].append((int(pos), row_str))

        for contig in contigs_order_list:
            row_list = [item[1] for item in sorted(contig_dict[contig], key=lambda x: x[0])]
            output.write(''.join(row_list))

    output.close()

    if genotyping_mode_vcf_fn is not None:
        print("[INFO] Total variants for genotyping: {}, total somatic variant calls: {}, added {} variants into output VCF"\
              .format(len(variant_dict), len(somatic_variant_dict), count))
    elif hybrid_mode_vcf_fn is not None:
        print("[INFO] Total additional variants for genotyping: {}, total somatic variant calls: {}, added {} variants into output VCF"\
              .format(len(variant_dict), len(somatic_variant_dict), count))

    compress_index_vcf(output_fn)


def main():
    parser = ArgumentParser(description="Genotype VCF in postprocessing")

    parser.add_argument('--genotyping_mode_vcf_fn', type=str_none, default=None,
                        help="Candidate sites VCF file input for genotyping")

    parser.add_argument('--hybrid_mode_vcf_fn', type=str_none, default=None,
                        help="Variants that passed the threshold and additional VCF candidates will both be subjected to variant calling")

    parser.add_argument('--candidates_folder', type=str, default=None,
                        help="Candidate folder to store the genotyped alternative candidate information")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="Somatic VCF input")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Output vcf file name")

    parser.add_argument('--switch_genotype', type=str2bool, default=True,
                        help="Switch missed variant genotype to ./.")

    args = parser.parse_args()

    genotype_vcf(args)


if __name__ == "__main__":
    main()
