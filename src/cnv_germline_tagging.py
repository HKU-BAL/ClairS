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
import sys
import os

from argparse import ArgumentParser
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__))))
from shared.vcf import VcfReader
from src.sort_vcf import compress_index_vcf
from numpy import *
from scipy.stats import binomtest

file_directory = os.path.dirname(os.path.realpath(__file__))
entry_path = os.path.join(file_directory, 'verdict')

def tumor_allele_counter_command(args):
    #allele counter
    command = f'time {args.parallel} -j{args.threads} '
    command += f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{args.allele_counter}/lib" {args.allele_counter}/bin/alleleCounter '
    command += f'-b {args.tumor_bam_fn} '
    command += f'-l {args.cnv_resource_dir}/loci_files/G1000_loci_hg38_{{1}}.txt '
    command += f'-o {args.output_dir}/{args.tumor_sample_name}_AlleleCount_{{1}}.txt '
    command += '-m 20 '
    command += '-q 20 '
    command += '-f 0 '
    command += '-F 2316 '
    command += '--dense-snps '
    command += f':::: {args.contig_fn}'

    return command

def normal_allele_counter_command(args):
    command = f'time {args.parallel} -j{args.threads} '
    command += f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{args.allele_counter}/lib" {args.allele_counter}/bin/alleleCounter '
    command += f'-b {args.normal_bam_fn} '
    command += f'-l {args.cnv_resource_dir}/loci_files/G1000_loci_hg38_{{1}}.txt '
    command += f'-o {args.output_dir}/{args.normal_sample_name}_AlleleCount_{{1}}.txt '
    command += '-m 20 '
    command += '-q 20 '
    command += '-f 0 '
    command += '-F 2316 '
    command += '--dense-snps '
    command += f':::: {args.contig_fn}'

    return command

def get_logr_baf_command(args):

    command = f'time {args.python} {args.verdict}/get_logr_and_baf.py '
    command += f'--tumor_allele_counts_file_prefix {args.output_dir}/{args.tumor_sample_name}_AlleleCount_ '
    command += f'--normal_allele_counts_file_prefix {args.output_dir}/{args.normal_sample_name}_AlleleCount_ '
    command += f'--alleles_file_prefix {args.cnv_resource_dir}/allele_files/G1000_alleles_hg38_ '
    command += f'--tumor_logr_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR.txt '
    command += f'--tumor_baf_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF.txt '
    command += f'--normal_baf_output_file {args.output_dir}/{args.normal_sample_name}_Normal_BAF.txt '
    command += f'--sample_name {args.tumor_sample_name} '
    command += f'--normal_sample_name {args.normal_sample_name}'

    return command

def correct_logr_command(args):
    command = f'time {args.python} {args.verdict}/correct_logr.py '
    command += f'--tumor_logr_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR.txt '
    command += f'--gc_content_file {args.cnv_resource_dir}/GC_G1000_hg38.txt '
    command += f'--replication_timing_file {args.cnv_resource_dir}/RT_G1000_hg38.txt '
    command += f'--tumor_logr_correction_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_Correction.txt '
    command += f'--sample_name {args.tumor_sample_name}'

    return command


def predict_germline_genotypes_command(args):
    command = f'time {args.python} {args.verdict}/predict_germline_genotypes.py '
    command += f'--tumor_logr_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_Correction.txt '
    command += f'--tumor_baf_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF.txt '
    command += f'--normal_baf_file {args.output_dir}/{args.normal_sample_name}_Normal_BAF.txt '
    command += f'--germline_genotypes_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_GG.txt '
    command += f'--sample_name {args.tumor_sample_name}'

    return command


def create_aspcf_command(args):
    command = f'time {args.python} {args.verdict}/aspcf.py '
    command += f'--tumor_logr_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_Correction.txt '
    command += f'--tumor_baf_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF.txt '
    command += f'--germline_genotypes_file {args.output_dir}/{args.tumor_sample_name}_Tumor_GG.txt '
    command += f'--tumor_logr_pcfed_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_PCFed.txt '
    command += f'--tumor_baf_pcfed_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF_PCFed.txt '
    command += f'--penalty 1000 '
    command += f'--sample_name {args.tumor_sample_name}'

    return command


def create_verdict_run_ascat_command(args):
    command = f'time {args.python} {args.verdict}/run_ascat.py '
    command += f'--tumor_logr_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_Correction.txt '
    command += f'--tumor_baf_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF.txt '
    command += f'--germline_genotypes_file {args.output_dir}/{args.tumor_sample_name}_Tumor_GG.txt '
    command += f'--tumor_logr_segmented_file {args.output_dir}/{args.tumor_sample_name}_Tumor_LogR_PCFed.txt '
    command += f'--tumor_baf_segmented_file {args.output_dir}/{args.tumor_sample_name}_Tumor_BAF_PCFed.txt '
    command += f'--tumor_purity_ploidy_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_Purity_Ploidy.txt '
    command += f'--tumor_cna_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_CNA.txt '
    command += f'--gamma 1.0 '
    command += f'--min_ploidy 1.5 '
    command += f'--max_ploidy 5.5 '
    command += f'--min_purity 0.1 '
    command += f'--max_purity 1.05 '
    command += f'--sample_name {args.tumor_sample_name}'

    return command

def tag_germline_variant(args):
    command = f'time {args.python} {args.verdict}/tag_germline_variant.py '
    command += f'--input_vcf_fn {args.input_vcf_fn} '
    command += f'--output_fn {args.output_fn} '
    command += f'--tumor_purity_ploidy_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_Purity_Ploidy.txt '
    command += f'--tumor_cna_output_file {args.output_dir}/{args.tumor_sample_name}_Tumor_CNA.txt '

    return command

def get_cnv_purity(args):

    args.verdict = args.verdict if args.verdict else entry_path

    if args.output_dir is not None and os.path.exists(args.output_dir):
        subprocess.run(f'mkdir -p {args.output_dir}', shell=True)

    tac_command = tumor_allele_counter_command(args)
    nac_command = normal_allele_counter_command(args)

    glb_command = get_logr_baf_command(args)
    cl_command = correct_logr_command(args)
    pgg_command = predict_germline_genotypes_command(args)
    ca_command = create_aspcf_command(args)
    cvra_command = create_verdict_run_ascat_command(args)
    tg_command = tag_germline_variant(args)

    commands_list = (tac_command, nac_command, glb_command, cl_command, pgg_command, ca_command, cvra_command, tg_command)
    for i, command in enumerate(zip(commands_list)):

        print(f"[INFO] STEP {i} RUN THE FOLLOWING COMMAND FOR CNV GERMLINE TAGGING:")
        print(command)
        try:
            return_code = subprocess.check_call(command, shell=True, stdout=sys.stdout)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("ERROR in STEP {}, THE FOLLOWING COMMAND FAILED: {}\n".format(i + 1, command))
            exit(1)
        print('')

def main():
    parser = ArgumentParser(description="add cnv tag into output_vcf file with the same input prefix")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input, required")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input directory")

    parser.add_argument('--tumor_purity', type=float, default=None,
                        help="Output directory")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output directory")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Output file name")

    parser.add_argument('--is_snv', action='store_true',
                        help="SNV input candidates")

    parser.add_argument('--is_indel', action='store_true',
                        help="Indel input_candidates")

    parser.add_argument('--python', type=str, default='python3',
                        help="Absolute path of python, python3 >= 3.9 is required")

    parser.add_argument('--parallel', type=str, default='parallel',
                        help="Absolute path of parallel, parallel >= 20191122 is required")

    parser.add_argument('--pypy', type=str, default='pypy',
                        help="Absolute path of pypy3, pypy3 >= 3.6 is required")

    parser.add_argument('--contig_fn', type=str, default=None,
                        help="All contig name")

    parser.add_argument('--tumor_sample_name', type=str, default='tumor',
                        help="Sample name")

    parser.add_argument('--normal_sample_name', type=str, default='normal',
                        help="Normal sample name")

    parser.add_argument('--allele_counter', type=str, default=None,
                        help="Directory of allele counter")

    parser.add_argument('--cnv_resource_dir', type=str, default=None,
                        help="Reference resource of CNV directory")

    parser.add_argument('--threads', type=int, default=None,
                        help="Threads to use")

    parser.add_argument('--verdict', type=str, default=None,
                        help="Path of entry path")

    args = parser.parse_args()

    get_cnv_purity(args)


if __name__ == "__main__":
    main()

