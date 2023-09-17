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
import os

from argparse import ArgumentParser


def concat_files(args):
    input_dir = args.input_dir
    input_prefix = args.input_prefix
    output_dir = args.output_dir if args.output_dir is not None else input_dir
    output_fn = args.output_fn
    is_snv = args.is_snv
    is_indel = args.is_indel

    if not os.path.exists(input_dir):
        sys.exit("[ERROR] The input prefix is not found: {}".format(input_prefix))

    if output_fn is not None and '/' not in output_fn:
        output_fn = os.path.join(output_dir, output_fn)

    if is_snv and output_fn is None:
        output_fn = os.path.join(output_dir, "CANDIDATES_FILES")
    elif is_indel and output_fn is None:
        output_fn = os.path.join(output_dir, "INDEL_CANDIDATES_FILES")

    with open(output_fn, 'w') as f:
        for file in os.listdir(input_dir):
            if file.startswith(input_prefix):
                tmp_f = open(os.path.join(input_dir, file))
                for row in tmp_f:
                    if row.rstrip() == '':
                        continue
                    f.write(row)
                tmp_f.close()



def main():
    parser = ArgumentParser(description="Concat file with the same input prefix")

    parser.add_argument('--input_dir', type=str, default=None, required=True,
                        help="Input directory")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output directory")

    parser.add_argument('--input_prefix', type=str, default=None, required=True,
                        help="Prefix of input file")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Output file name")

    parser.add_argument('--is_snv', action='store_true',
                        help="SNV input candidates")

    parser.add_argument('--is_indel', action='store_true',
                        help="Indel input_candidates")

    args = parser.parse_args()

    concat_files(args)

if __name__ == "__main__":
    main()
