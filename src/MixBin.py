import os
import random
import shlex
import subprocess
from argparse import ArgumentParser, SUPPRESS
from subprocess import PIPE, run as subprocess_run
random.seed(0)

from src.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool

cov_suffix = ".depth.mosdepth.summary.txt"

def get_coverage(depth_log):
    last_row = open(depth_log).readlines()[-1]
    depth = int(float(last_row.split()[3]))
    return depth

def MixBin(args):
    normal_bam_fn = args.normal_bam_fn
    # tumor_bam_fn = args.tumor_bam_fn
    output_fn = args.output_fn
    input_dir = args.input_dir
    cov_dir = args.cov_dir
    ctg_name = args.ctgName
    samtools_execute_command = args.samtools
    samtools_threads = args.samtools_threads
    min_coverage = args.minCoverage
    synthetic_proportion = args.synthetic_proportion
    contaminative_proportion = args.contaminative_proportion
    platform = args.platform
    normal_bam_depth = args.normal_bam_depth
    tumor_bam_depth = args.tumor_bam_depth

    normal_depth_log = os.path.join(cov_dir, 'normal_' + ctg_name + cov_suffix)
    tumor_depth_log = os.path.join(cov_dir, 'tumor_' + ctg_name + cov_suffix)
    normal_bam_depth = get_coverage(normal_depth_log)
    tumor_bam_depth = get_coverage(tumor_depth_log)

    normal_bin_num = int(int(normal_bam_depth) / int(min_coverage))
    tumor_bin_num = int(int(tumor_bam_depth) / int(min_coverage))

    bam_list = os.listdir(input_dir)
    normal_bam_list = [bam for bam in bam_list if bam.startswith('normal_' + ctg_name + '_')]
    tumor_bam_list = [bam for bam in bam_list if bam.startswith('tumor_' + ctg_name  + '_')]

    tumor_depth = (normal_bam_depth) * synthetic_proportion
    normal_depth = normal_bam_depth - tumor_depth
    sampled_normal_bin_num = int(normal_depth / min_coverage)
    sampled_tumor_bin_num = int(tumor_depth / min_coverage)

    random.seed(0)
    sampled_normal_bam_list = [normal_bam_list[idx] for idx in random.sample(range(normal_bin_num), sampled_normal_bin_num)]
    sampled_tumor_bam_list = [tumor_bam_list[idx] for idx in random.sample(range(tumor_bin_num), sampled_tumor_bin_num)]

    input_sampled_bam_list = [bam for bam in sampled_normal_bam_list + sampled_tumor_bam_list]
    print (input_sampled_bam_list)
    input_sampled_bam_list = ' '.join([os.path.join(input_dir, bam) for bam in input_sampled_bam_list])

    subprocess_run(shlex.split("{} merge -f -@{} {} {}".format(samtools_execute_command, samtools_threads,output_fn, input_sampled_bam_list)))
    subprocess_run(shlex.split("{} index -@{} {}".format(samtools_execute_command, samtools_threads,output_fn)))

    if contaminative_proportion is None:
        normal_output_bam = output_fn.replace('tumor_', 'normal_')
        subprocess_run(shlex.split("ln -sf {} {}".format(normal_bam_fn, normal_output_bam)))
        subprocess_run(shlex.split("ln -sf {}.bai {}.bai".format(normal_bam_fn, normal_output_bam)))

    #
    # for bam_fn, bin_num, prefix in zip((normal_bam_fn, tumor_bam_fn), (normal_bin_num, tumor_bin_num), ("normal", 'tumor')):
    #     subprocess_list = []
    #     for bin_idx in range(bin_num):
    #         output_fn = os.path.join(output_dir, prefix + "_" + str(bin_idx))
    #         save_file_fp = subprocess_popen(
    #             shlex.split("{} view -bh - -o {}".format(samtools_execute_command, output_fn)), stdin=PIPE,
    #             stdout=PIPE)
    #         subprocess_list.append(save_file_fp)
    #
    #     samtools_view_command = "{} view -@ {} -h {} {}".format(samtools_execute_command, samtools_threads, bam_fn, ctg_name if ctg_name else "")
    #
    #     samtools_view_process = subprocess_popen(shlex.split(samtools_view_command))
    #
    #     for row_id, row in enumerate(samtools_view_process.stdout):
    #         if row[0] == '@':
    #             for subprocess in subprocess_list:
    #                 subprocess.stdin.write(row)
    #             continue
    #         bin_id = int(random.random() * 100) % bin_num
    #         subprocess_list[bin_id].stdin.write(row)
    #
    #     samtools_view_process.stdout.close()
    #     samtools_view_process.wait()
    #     for save_file_fp in subprocess_list:
    #         save_file_fp.stdin.close()
    #         save_file_fp.wait()


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--normal_bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--tumor_bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--normal_bam_depth', type=int, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--tumor_bam_depth', type=int, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--input_dir', type=str, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--minCoverage', type=int, default=4,
                        help="Reference fasta file input, required")

    parser.add_argument('--samtools_threads', type=int, default=32,
                        help="Reference fasta file input, required")

    parser.add_argument('--synthetic_proportion', type=float, default=0.25,
                        help="Reference fasta file input, required")

    parser.add_argument('--contaminative_proportion', type=float, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--cov_dir', type=str, default=None,
                        help="Sorted BAM file input, required")


    args = parser.parse_args()

    MixBin(args)


if __name__ == "__main__":
    main()
