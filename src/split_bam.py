import os
import random
import shlex

from argparse import ArgumentParser, SUPPRESS
from subprocess import PIPE

from src.utils import subprocess_popen

random.seed(0)
cov_suffix = ".cov.mosdepth.summary.txt"


def get_coverage(coverage_log, ctg_name=None):
    # we use the overall average coverage if no contig specified
    if ctg_name is None:
        last_row = open(coverage_log).readlines()[-1]
        coverage = int(float(last_row.split()[3]))
    else:
        all_rows = open(coverage_log).readlines()
        ctg_row = [row for row in all_rows if row.split()[0] == ctg_name]
        if len(ctg_row) == 0:
            print('[ERROR] no contig coverage found for contig {}'.format(ctg_name))
        coverage = int(ctg_row[0].split()[3])
    return coverage


def split_bam(args):
    normal_bam_fn = args.normal_bam_fn
    tumor_bam_fn = args.tumor_bam_fn
    output_dir = args.output_dir
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    samtools_threads = args.samtools_threads
    samtools_output_threads = args.samtools_output_threads
    min_bin_coverage = args.min_bin_coverage
    cov_dir = args.cov_dir

    normal_coverage_log = os.path.join(cov_dir, 'raw_normal_' + ctg_name + cov_suffix)
    tumor_coverage_log = os.path.join(cov_dir, 'raw_tumor_' + ctg_name + cov_suffix)
    normal_bam_coverage = get_coverage(normal_coverage_log)
    tumor_bam_coverage = get_coverage(tumor_coverage_log)

    normal_bin_num = int(int(normal_bam_coverage) / int(min_bin_coverage))
    tumor_bin_num = int(int(tumor_bam_coverage) / int(min_bin_coverage))

    for bam_fn, bin_num, prefix in zip((normal_bam_fn, tumor_bam_fn), (normal_bin_num, tumor_bin_num),
                                       ("normal", 'tumor')):
        subprocess_list = []
        for bin_idx in range(bin_num):
            output_fn = os.path.join(output_dir, '_'.join([prefix, ctg_name, str(bin_idx)]) + '.bam')
            save_file_fp = subprocess_popen(
                shlex.split(
                    "{} view -bh -@ {} - -o {}".format(samtools_execute_command, samtools_output_threads, output_fn)),
                stdin=PIPE,
                stdout=PIPE)
            subprocess_list.append(save_file_fp)

        samtools_view_command = "{} view -@ {} -h {} {}".format(samtools_execute_command, samtools_threads, bam_fn,
                                                                ctg_name if ctg_name else "")

        samtools_view_process = subprocess_popen(shlex.split(samtools_view_command))

        for row_id, row in enumerate(samtools_view_process.stdout):
            if row[0] == '@':
                for subprocess in subprocess_list:
                    subprocess.stdin.write(row)
                continue
            # replace random shuffle to partition for reproducibility
            # bin_id = int(random.random() * 100) % bin_num
            bin_id = row_id % bin_num
            # add prefix  for each normal and tumor reads
            subprocess_list[bin_id].stdin.write(prefix[0] + row)

        samtools_view_process.stdout.close()
        samtools_view_process.wait()
        for save_file_fp in subprocess_list:
            save_file_fp.stdin.close()
            save_file_fp.wait()

    print("[INFO] Contig/Normal coverage/Tumor coverage: {}/{}/{}".format(ctg_name, normal_bam_coverage,
                                                                          tumor_bam_coverage))


def main():
    parser = ArgumentParser(description="Split BAM into small bins")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Sorted chunked BAM file output directory, required")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input, required")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--normal_bam_coverage', type=int, default=None,
                        help="Normal BAM coverage calculated using mosdepth")

    parser.add_argument('--tumor_bam_coverage', type=int, default=None,
                        help="Tumor BAM coverage calculated using mosdepth")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required")

    parser.add_argument('--cov_dir', type=str, default=None,
                        help="Directory of mosdepth coverage summary")

    parser.add_argument('--min_bin_coverage', type=int, default=4,
                        help="Minimum average coverage for each chunked bin")

    parser.add_argument('--samtools_threads', type=int, default=32,
                        help="Samtools threads to read input BAM")

    parser.add_argument('--samtools_output_threads', type=int, default=24,
                        help="Samtools threads to write input BAM")

    args = parser.parse_args()

    split_bam(args)


if __name__ == "__main__":
    main()
