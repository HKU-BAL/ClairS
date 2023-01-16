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

import sys
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
import os

from shared.utils import subprocess_popen


class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def remove_common_suffix(ref_base, alt_base):
    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))  # keep at least one base
    prefix = ref_base[::-1]
    for string in alt_base:
        string = string[::-1]
        while string[:len(prefix)] != prefix and prefix:
            prefix = prefix[:len(prefix) - 1]
        if not prefix:
            break
    res_length = len(prefix)
    if res_length > min_length:
        return ref_base, alt_base
    return ref_base[:len(ref_base) - res_length], [item[:len(item) - res_length] for item in alt_base]

    return ref_base[-min_length], [item[-min_length] for item in alt_base]


def decode_alt(ref_base, alt_base):
    if ',' not in alt_base:
        return [ref_base], [alt_base]
    alt_base = alt_base.split(',')
    ref_base_list, alt_base_list = [], []
    for ab in alt_base:
        rb,ab = remove_common_suffix(ref_base, [ab])
        ref_base_list.append(rb)
        alt_base_list.append(ab[0])
    return ref_base_list, alt_base_list

def UpdateVar(args):
    var_fn = args.var_fn
    input_var_fn = args.input_var_fn
    ctg_name = args.ctgName
    alt_fn_prefix = args.alt_fn_prefix

    if args.var_fn != "PIPE":
        var_fpo = open(var_fn, "wb")
        var_fp = subprocess_popen(shlex.split("gzip -c"), stdin=PIPE, stdout=var_fpo)
    else:
        var_fp = TruthStdout(sys.stdout)

    input_var_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (input_var_fn)))

    Y = {}
    for row in input_var_fp.stdout:
        if row[0] == "#":
            continue
        columns = row.strip().split()
        contig_name, pos, ref_base, alt_base, genotype1, genotype2 = columns
        if ctg_name and contig_name != ctg_name:
            continue
        ref_base_list, alt_base_list = decode_alt(ref_base, alt_base)
        Y[pos] = [ref_base_list, alt_base_list, row]

    alt_fn_prefix = alt_fn_prefix.split('/')
    directry, file_prefix = '/'.join(alt_fn_prefix[:-1]), alt_fn_prefix[-1]
    file_list = [f for f in os.listdir(directry) if f.startswith(file_prefix)]

    found_num = 0
    for f in file_list:
        for row_id, row in enumerate(open(os.path.join(directry, f), 'r')):
            # print (row_id)
            chr_pos, depth, alt_info = row.split('\t')[:3]
            contig_name, pos = chr_pos.split()
            if ctg_name and contig_name != ctg_name:
                continue
            if pos not in Y:
                continue
            seqs = alt_info.split(' ')
            seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
            # alt_list = sorted(list(seq_alt_bases_dict.items()), key=lambda x: x[1], reverse=True)

            ref_base_list, alt_base_list, r = Y[pos]
            found = 0
            for alt_type in seq_alt_bases_dict:
                if '*' in alt_type or '#' in alt_type or 'R' in alt_type:
                    continue
                if alt_type[0] == 'X':
                    if alt_type[1] in alt_base_list:
                        found += 1
                elif alt_type[0] == 'I':
                    if alt_type[1:] in alt_base_list:
                        found += 1
                elif alt_type[0] == 'D':
                    del_cigar = alt_type[1:]
                    for rb, ab in zip(ref_base_list, alt_base_list):
                        if rb[1:] == del_cigar and len(ab) == 1:
                            found += 1
            if found >= len(alt_base_list):
                found_num +=1
                var_fp.stdin.write(r)
            else:
                # print (r)
                var_fp.stdin.write(" ".join(r.split(' ')[:4]) +' -1 -1\n')
    print("total variants/find variants/ratio: {}/{}/{}".format(len(Y), found_num,
                                                                   round(found_num/float(len(Y)), 4)),
          file=sys.stderr)
    input_var_fp.stdout.close()
    input_var_fp.wait()

    if args.var_fn != "PIPE":
        var_fp.stdin.close()
        var_fp.wait()
        var_fpo.close()

def main():
    parser = ArgumentParser(description="Extract variant type and allele from a Truth dataset")

    parser.add_argument('--vcf_fn', type=str, default="input.vcf",
                        help="Truth vcf file input, default: %(default)s")

    parser.add_argument('--var_fn', type=str, default="PIPE",
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--alt_fn_prefix', type=str, default=None,
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--input_var_fn', type=str, default="PIPE",
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    UpdateVar(args)


if __name__ == "__main__":
    main()
