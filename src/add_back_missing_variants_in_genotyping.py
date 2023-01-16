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

from sys import stderr
from subprocess import PIPE
from argparse import ArgumentParser
from collections import defaultdict
from subprocess import Popen

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23))] + [str(a) for a in
                                                                                   list(range(1, 23))]


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


def genotype_vcf(args):
    vcf_fn = args.vcf_fn
    call_fn = args.call_fn
    output_fn = args.output_fn
    switch_genotype = args.switch_genotype

    output = open(output_fn, 'w')
    vcf_reader = VcfReader(vcf_fn=vcf_fn,
                           ctg_name=None,
                           keep_row_str=True,
                           save_header=True)

    vcf_reader.read_vcf()
    variant_dict = vcf_reader.variant_dict

    rc = subprocess.run('cp {} {}.bak'.format(call_fn, call_fn), shell=True)

    somatic_vcf_reader = VcfReader(vcf_fn=call_fn,
                                  ctg_name=None,
                                  keep_row_str=True,
                                  save_header=True)

    somatic_vcf_reader.read_vcf()
    somatic_variant_dict = somatic_vcf_reader.variant_dict

    all_contigs_list = list(set([k[0] for k in variant_dict]))

    contigs_order = major_contigs_order + all_contigs_list

    contigs_order_list = sorted(all_contigs_list, key=lambda x: contigs_order.index(x))

    output.write(somatic_vcf_reader.header)



    count = 0
    contig_dict = defaultdict(list)
    for k, v in variant_dict.items():
        ctg, pos = k
        if k not in somatic_variant_dict:
            row_str = variant_dict[k].row_str
            count += 1
            if switch_genotype:
                columns = row_str.rstrip().split('\t')
                if len(columns) < 10:
                    columns += ['.'] * (10 - len(columns))
                columns[3] = columns[3][0] if len(columns[3]) > 0 else '.' # Only keep the reference base for REF
                columns[4] = '.' # ALT to 0
                columns[5] = "." # QUAL to .
                columns[6] = "." # FILTER to .
                columns[7] = "." # INFO to .
                columns[8] = "GT" # keep GT tag only
                columns[9] = './.'
                row_str = '\t'.join(columns) + '\n'
        else:
            row_str = somatic_variant_dict[k].row_str

        contig_dict[ctg].append((int(pos), row_str))

    for contig in contigs_order_list:
        row_list = [item[1] for item in sorted(contig_dict[contig], key=lambda x: x[0])]
        output.write(''.join(row_list))

    output.close()

    print("[INFO] Total variants for genotyping: {}, total somatic variant calls: {}, added {} variants into output VCF"\
          .format(len(variant_dict), len(somatic_variant_dict), count))
    compress_index_vcf(output_fn)


def main():
    parser = ArgumentParser(description="Genotype VCF in postprocessing")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input for genotyping")

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
