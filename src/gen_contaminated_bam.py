import os
import random
import shlex
import subprocess

from argparse import ArgumentParser, SUPPRESS
from subprocess import run as subprocess_run

from src.utils import str2bool, str_none

random.seed(0)
cov_suffix = ".cov.mosdepth.summary.txt"

def get_coverage_from_bam(args, bam_fn, is_tumor=False, output_prefix=None):
    ctg_name = args.ctg_name
    dry_run = args.dry_run

    mosdepth = args.mosdepth
    contig_option = "" if ctg_name is None else "-c {}".format(ctg_name)

    prefix = 'tumor' if is_tumor else 'normal'
    if output_prefix is None:
        output_prefix = os.path.join(args.output_dir, 'tmp', 'raw_{}'.format(prefix))


    mos_depth_command = "{} -t {} {} -n -x --quantize 0:15:150: {}.cov {}".format(mosdepth,
                                                                 args.samtools_threads,
                                                                 contig_option,
                                                                 output_prefix,
                                                                 bam_fn)

    print("[INFO] Calculating coverge for {} BAM using mosdepth...".format(prefix))


    print('[INFO] Will run the following commands:')
    print(mos_depth_command)
    if not dry_run:
        subprocess.run(mos_depth_command, shell=True)

    coverage_log = os.path.join(output_prefix + cov_suffix)
    if ctg_name is None:
        last_row = open(coverage_log).readlines()[-1]
        coverage = float(float(last_row.split()[3]))
    else:
        all_rows = open(coverage_log).readlines()
        ctg_row = [row for row in all_rows if row.split()[0] == ctg_name]
        if len(ctg_row) == 0:
            print('[ERROR] no contig coverage found for contig {}'.format(ctg_name))
        coverage = float(ctg_row[0].split()[3])
    return coverage


def random_sample(population, k, seed=0):
    random.seed(seed)
    return random.sample(population, k)


def gen_contaminated_bam(args):
    tumor_bam_fn = args.tumor_bam_fn
    normal_bam_fn = args.normal_bam_fn
    output_dir = args.output_dir
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    samtools_threads = args.samtools_threads
    normal_purity = args.normal_purity
    tumor_purity = args.tumor_purity

    if not os.path.exists(output_dir):
        rc = subprocess.run('mkdir -p {}/tmp'.format(output_dir), shell=True)
    normal_bam_coverage = args.normal_bam_coverage if args.normal_bam_coverage else get_coverage_from_bam(args,normal_bam_fn, False)
    tumor_bam_coverage = args.tumor_bam_coverage if args.tumor_bam_coverage else get_coverage_from_bam(args, tumor_bam_fn, True)

    print("[INFO] Normal/Tumor BAM coverage: {}/{}".format(normal_bam_coverage,
                                                           tumor_bam_coverage))

    # add normal to normal
    if normal_purity is not None:
        contam_coverage = normal_bam_coverage * (1 - normal_purity)
        rest_normal_coverage = normal_bam_coverage - contam_coverage

        tumor_subsample_pro = "%.3f" % (contam_coverage / float(tumor_bam_coverage))
        normal_subsample_pro = "%.3f" % (rest_normal_coverage / float(normal_bam_coverage))

        print("[INFO] Normal/Tumor subsample proportion: {}/{}".format(normal_subsample_pro,
                                                               tumor_subsample_pro))

        tumor_subsample_bam = os.path.join(args.output_dir, 'tmp', 'tumor_subsample.bam')
        normal_subsample_bam = os.path.join(args.output_dir, 'tmp', 'normal_rest.bam')

        contig_option = "" if ctg_name is None else ctg_name

        #tumor subsample cmd
        t_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                          samtools_threads,
                                                          tumor_subsample_pro,
                                                          tumor_subsample_bam,
                                                          tumor_bam_fn,
                                                          contig_option,
                                                          )
        #normal subsample cmd
        n_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                          samtools_threads,
                                                          normal_subsample_pro,
                                                          normal_subsample_bam,
                                                          normal_bam_fn,
                                                          contig_option)

        n_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, tumor_subsample_bam)
        t_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, normal_subsample_bam)


        print('[INFO] Will run the following commands:')
        print(t_s_cmd)
        print(n_s_cmd)
        print(n_index_cmd)
        print(t_index_cmd)

        if not args.dry_run:
            subprocess.run(t_s_cmd, shell=True)
            subprocess.run(n_s_cmd, shell=True)
            subprocess.run(n_index_cmd, shell=True)
            subprocess.run(t_index_cmd, shell=True)

        normal_output_bam = os.path.join(output_dir, "normal_purity_{}.bam".format(normal_purity))

        print("[INFO] Merging tumor BAM into normal BAM as contamination...")
        merge_cmd = "{} merge -f -@{} {} {} {}".format(samtools_execute_command, samtools_threads, normal_output_bam,
                                            tumor_subsample_bam, normal_subsample_bam)
        index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, normal_output_bam)


        print('[INFO] Will run the following commands:')
        print(merge_cmd)
        print(index_cmd)
        if not args.dry_run:
            subprocess_run(merge_cmd, shell=True)
            subprocess.run(index_cmd, shell=True)

        if args.cal_output_bam_coverage:
            get_coverage_from_bam(args, normal_output_bam, False, os.path.join(args.output_dir, 'cov'))

        if args.remove_intermediate_dir:
            tmp_file_path = os.path.join(args.output_dir, 'tmp')
            if os.path.exist(tmp_file_path):
                rc = subprocess.run("rm -rf {}".format(tmp_file_path), shell=True)

        print("[INFO] Finishing merging, output file: {}".format(normal_output_bam))

    #add normal to tumor
    if tumor_purity is not None:
        contam_coverage = tumor_bam_coverage * (1 - tumor_purity)
        rest_tumor_coverage = tumor_bam_coverage - contam_coverage

        tumor_subsample_pro = "%.3f" % (rest_tumor_coverage / float(tumor_bam_coverage))
        normal_subsample_pro = "%.3f" % (contam_coverage / float(normal_bam_coverage))

        print("[INFO] Normal/Tumor subsample proportion: {}/{}".format(normal_subsample_pro,
                                                               tumor_subsample_pro))

        tumor_subsample_bam = os.path.join(args.output_dir, 'tmp', 'tumor_rest.bam')
        normal_subsample_bam = os.path.join(args.output_dir, 'tmp', 'normal_subsample.bam')

        contig_option = "" if ctg_name is None else ctg_name

        #tumor subsample cmd
        t_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                          samtools_threads,
                                                          tumor_subsample_pro,
                                                          tumor_subsample_bam,
                                                          tumor_bam_fn,
                                                          contig_option,
                                                          )
        #normal subsample cmd
        n_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                          samtools_threads,
                                                          normal_subsample_pro,
                                                          normal_subsample_bam,
                                                          normal_bam_fn,
                                                          contig_option)

        n_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, tumor_subsample_bam)
        t_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, normal_subsample_bam)


        print('[INFO] Will run the following commands:')
        print(t_s_cmd)
        print(n_s_cmd)
        print(n_index_cmd)
        print(t_index_cmd)

        if not args.dry_run:
            subprocess.run(t_s_cmd, shell=True)
            subprocess.run(n_s_cmd, shell=True)
            subprocess.run(n_index_cmd, shell=True)
            subprocess.run(t_index_cmd, shell=True)

        tumor_output_bam = os.path.join(output_dir, "tumor_contaminated_{}.bam".format(tumor_purity))

        print("[INFO] Merging normal BAM into tumor BAM as contamination...")
        merge_cmd = "{} merge -f -@{} {} {} {}".format(samtools_execute_command, samtools_threads, tumor_output_bam,
                                            tumor_subsample_bam, normal_subsample_bam)
        index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, tumor_output_bam)

        print('[INFO] Will run the following commands:')
        print(merge_cmd)
        print(index_cmd)
        if not args.dry_run:
            subprocess_run(merge_cmd, shell=True)
            subprocess.run(index_cmd, shell=True)

        if args.cal_output_bam_coverage:
            get_coverage_from_bam(args, tumor_output_bam, True, os.path.join(args.output_path, 'cov'))

        if args.remove_intermediate_dir:
            tmp_file_path = os.path.join(args.output_dir, 'tmp')
            if os.path.exist(tmp_file_path):
                rc = subprocess.run("rm -rf {}".format(tmp_file_path), shell=True)

        print("[INFO] Finishing merging, output file: {}".format(tumor_output_bam))


def main():
    parser = ArgumentParser(description="Generate contaminated BAM for benchmarking")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--normal_bam_coverage', type=int, default=None,
                        help="Normal BAM coverage calculated using mosdepth")

    parser.add_argument('--tumor_bam_coverage', type=int, default=None,
                        help="Tumor BAM coverage calculated using mosdepth")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Sorted chunked BAM file output path")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required")

    parser.add_argument('--samtools_threads', type=int, default=32,
                        help="Samtools threads to read input BAM")

    parser.add_argument('--normal_purity', type=float, default=None,
                        help="Normal purity")

    parser.add_argument('--tumor_purity', type=float, default=None,
                        help="Tumor purity")

    parser.add_argument('--ctg_name', type=str_none, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--mosdepth', type=str, default="mosdepth",
                        help="Path to the 'mosdepth'")

    parser.add_argument('--dry_run', type=str2bool, default=0,
                        help="EXPERIMENTAL: Only print the synthetic log, debug only")

    parser.add_argument('--cal_output_bam_coverage', type=str2bool, default=0,
                        help="EXPERIMENTAL: Only print the synthetic log, debug only")

    parser.add_argument('--remove_intermediate_dir', action='store_true',
                        help="Remove intermediate directory. Default: False")

    args = parser.parse_args()

    gen_contaminated_bam(args)


if __name__ == "__main__":
    main()
