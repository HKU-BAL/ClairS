import os
import random
import shlex

from argparse import ArgumentParser, SUPPRESS
from subprocess import run as subprocess_run

from src.utils import str2bool

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


def check_max_sampled_coverage(nor_cov, tum_cov, synthetic_proportion, pair_gamma=0.5, min_bin_coverage=4):
    # nor_cov = int(nor_cov / min_bin_coverage * min_bin_coverage)
    # tum_cov = int(tum_cov / min_bin_coverage * min_bin_coverage)
    max_synthetic_coverage_for_tumor = int(tum_cov / synthetic_proportion)
    max_synthetic_coverage_for_normal = int(nor_cov / (1 + pair_gamma - synthetic_proportion))

    max_synthetic_coverage = min(max_synthetic_coverage_for_tumor, max_synthetic_coverage_for_normal)
    max_synthetic_coverage = int(max_synthetic_coverage)

    return max_synthetic_coverage


def random_sample(population, k, seed=0):
    random.seed(seed)
    return random.sample(population, k)


def MixBin(args):
    tumor_bam_fn = args.tumor_bam_fn
    output_fn = args.output_fn
    input_dir = args.input_dir
    cov_dir = args.cov_dir
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    samtools_threads = args.samtools_threads
    tensor_sample_mode = args.tensor_sample_mode
    min_bin_coverage = args.min_bin_coverage
    synthetic_proportion = args.synthetic_proportion
    synthetic_coverage = args.synthetic_coverage
    contaminative_proportions = args.contaminative_proportions
    normal_coverage_proportion = args.normal_coverage_proportion
    normal_coverage_log = os.path.join(cov_dir, 'raw_normal_' + ctg_name + cov_suffix)
    tumor_coverage_log = os.path.join(cov_dir, 'raw_tumor_' + ctg_name + cov_suffix)
    normal_bam_coverage = args.normal_bam_coverage if args.normal_bam_coverage else get_coverage(normal_coverage_log)
    tumor_bam_coverage = args.tumor_bam_coverage if args.tumor_bam_coverage else get_coverage(tumor_coverage_log)

    normal_bin_num = int(int(normal_bam_coverage) / int(min_bin_coverage))
    tumor_bin_num = int(int(tumor_bam_coverage) / int(min_bin_coverage))

    if normal_coverage_proportion is not None and ctg_name is not None and type(int(ctg_name[3:])) == int:
        if int(ctg_name[3:]) <= 8:
            normal_coverage_proportion = 1.0
        elif int(ctg_name[3:]) > 16:
            normal_coverage_proportion = 0.5
        else:
            normal_coverage_proportion = 0.75

    bam_list = os.listdir(input_dir)
    if args.dry_run:
        normal_bam_list = ['normal_' + str(idx) for idx in range(normal_bin_num)]
        tumor_bam_list = ['tumor_' + str(idx) for idx in range(tumor_bin_num)]
    else:
        normal_bam_list = [bam for bam in bam_list if bam.startswith('normal_' + ctg_name + '_')]
        tumor_bam_list = [bam for bam in bam_list if bam.startswith('tumor_' + ctg_name + '_')]
    assert len(normal_bam_list) == normal_bin_num
    assert len(tumor_bam_list) == tumor_bin_num

    max_syn_coverage = check_max_sampled_coverage(nor_cov=normal_bam_coverage,
                                                  tum_cov=tumor_bam_coverage,
                                                  synthetic_proportion=synthetic_proportion,
                                                  pair_gamma=normal_coverage_proportion,
                                                  min_bin_coverage=min_bin_coverage)

    if synthetic_coverage is None:
        print('[INFO] Synthetic coverage not set, use maximum synthetic coverage {}'.format(max_syn_coverage))
        synthetic_coverage = max_syn_coverage
    elif synthetic_coverage > max_syn_coverage:
        print('[WARNING] Synthetic coverage is larger than maximum coverage {} > {}'.format(synthetic_coverage,
                                                                                            max_syn_coverage))
        synthetic_coverage = max_syn_coverage

    tumor_coverage = int(synthetic_coverage * synthetic_proportion)
    normal_coverage = int(synthetic_coverage - tumor_coverage)
    sampled_normal_bin_num = int(normal_coverage / min_bin_coverage)
    sampled_tumor_bin_num = int(tumor_coverage / min_bin_coverage)

    sampled_normal_bam_list = [normal_bam_list[idx] for idx in
                               random_sample(range(normal_bin_num), sampled_normal_bin_num,
                                             seed=int(synthetic_proportion * 100))]
    if tensor_sample_mode:
        sampled_tumor_bam_list = tumor_bam_list
    else:
        sampled_tumor_bam_list = [tumor_bam_list[idx] for idx in
                                  random_sample(range(tumor_bin_num), sampled_tumor_bin_num,
                                                seed=int(synthetic_proportion * 100))]

    # the rest sampled sample to generate normal bam
    rest_normal_bam_list = [normal_bam_fn for normal_bam_fn in normal_bam_list if
                            normal_bam_fn not in sampled_normal_bam_list]
    tumor_all_bin_num = sampled_normal_bin_num + sampled_tumor_bin_num
    normal_all_bin_num = int(tumor_all_bin_num * normal_coverage_proportion)
    if normal_all_bin_num <= len(rest_normal_bam_list):
        pair_normal_bam_list = [rest_normal_bam for rest_normal_bam in
                                random_sample(rest_normal_bam_list, normal_all_bin_num,
                                              seed=int(synthetic_proportion * 100))]
    else:
        # sample bin exhaust if pair normal coverage could not satisfy requirement
        # avoid to resample bin for model robustness
        print('[WARNING] Need to resample normal chunked BAMs!')
        resampled_bin_count = normal_all_bin_num - len(rest_normal_bam_list)
        resample_bin_list = [sampled_normal_bam_list[idx] for idx in
                             random_sample(range(len(sampled_normal_bam_list)), resampled_bin_count,
                                           seed=int(synthetic_proportion * 100))]
        pair_normal_bam_list = resample_bin_list + rest_normal_bam_list

    print(
        "[INFO] Raw normal BAM coverage/Raw tumor BAM coverage: {}x/{}x, normal sampled bins/tumor sampled bins:{}/{}".format(
            normal_bam_coverage, tumor_bam_coverage, sampled_normal_bin_num, sampled_tumor_bin_num))
    print("[INFO] Tumor sampled normal chunked BAMs coverage/bins: {}x/{}:{}".format(
        len(sampled_normal_bam_list) * min_coverage, len(sampled_normal_bam_list), ' '.join(sampled_normal_bam_list)))
    print("[INFO] Tumor sampled tumor chunked BAMs coverage/bins: {}x/{}:{}".format(
        len(sampled_tumor_bam_list) * min_coverage, len(sampled_tumor_bam_list), ' '.join(sampled_tumor_bam_list)))
    print("[INFO] Normal sampled BAMs coverage/bins: {}x/{}:{}".format(len(pair_normal_bam_list) * min_coverage,
                                                                       len(pair_normal_bam_list),
                                                                       ' '.join(pair_normal_bam_list)))
    print("[INFO] Synthetic coverage: {}, Normal sampled BAMs intersection: {}\n".format(synthetic_coverage, set(
        pair_normal_bam_list).intersection(set(sampled_normal_bam_list + sampled_tumor_bam_list))))

    print(tumor_bam_fn is not None and synthetic_proportion == 1.0 and len(sampled_normal_bam_list + sampled_tumor_bam_list) == tumor_bin_num \
            and len(sampled_normal_bam_list) == 0)

    if args.dry_run:
        return

    tumor_sampled_bam_list = sampled_normal_bam_list + sampled_tumor_bam_list
    tumor_sampled_bam_list = ' '.join([os.path.join(input_dir, bam) for bam in tumor_sampled_bam_list])
    tumor_output_bam = output_fn

    normal_output_bam = output_fn.replace('tumor_', 'normal_')
    normal_sampled_bam_list = ' '.join([os.path.join(input_dir, bam) for bam in pair_normal_bam_list])

    # no need to merge bins if use all tumor bins
    if tumor_bam_fn is not None and synthetic_proportion == 1.0 and len(tumor_sampled_bam_list) == tumor_bin_num \
            and len(sampled_normal_bam_list) == 0:
        subprocess_run(shlex.split("ln -sf {} {}".format(tumor_bam_fn, tumor_output_bam)))
        subprocess_run(shlex.split("ln -sf {}.bai {}.bai".format(tumor_bam_fn, tumor_output_bam)))

    else:
        subprocess_run(shlex.split(
            "{} merge -f -@{} {} {}".format(samtools_execute_command, samtools_threads, tumor_output_bam,
                                            tumor_sampled_bam_list)))
        subprocess_run(shlex.split("{} index -@{} {}".format(samtools_execute_command, samtools_threads, tumor_output_bam)))

    subprocess_run(shlex.split(
        "{} merge -f -@{} {} {}".format(samtools_execute_command, samtools_threads, normal_output_bam,
                                        normal_sampled_bam_list)))
    subprocess_run(
        shlex.split("{} index -@{} {}".format(samtools_execute_command, samtools_threads, normal_output_bam)))

    # if contaminative_proportions is not None:
    #     # contam_tumor_num =
    #     contaminative_proportion_list = contaminative_proportions.split(',')
    #     normal_output_bam = output_fn.replace('tumor_', 'normal_')
    #     normal_sampled_bam_list = ' '.join([os.path.join(input_dir, bam) for bam in pair_normal_bam_list])
    #     subprocess_run(shlex.split(
    #         "{} merge -f -@{} {} {}".format(samtools_execute_command, samtools_threads, normal_output_bam,
    #                                         normal_sampled_bam_list)))
    #     subprocess_run(
    #         shlex.split("{} index -@{} {}".format(samtools_execute_command, samtools_threads, normal_output_bam)))

    # subprocess_run(shlex.split("ln -sf {} {}".format(normal_bam_fn, normal_output_bam)))
    # subprocess_run(shlex.split("ln -sf {}.bai {}.bai".format(normal_bam_fn, normal_output_bam)))
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
    parser = ArgumentParser(description="Mix splited chunk BAM into multiple BAM")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--normal_bam_coverage', type=int, default=None,
                        help="Normal BAM coverage calculated using mosdepth")

    parser.add_argument('--tumor_bam_coverage', type=int, default=None,
                        help="Tumor BAM coverage calculated using mosdepth")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Sorted chunked BAM file output path")

    parser.add_argument('--input_dir', type=str, default=None,
                        help="Input directory, required")

    parser.add_argument('--min_bin_coverage', type=int, default=4,
                        help="Reference fasta file input, required")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required")

    parser.add_argument('--samtools_threads', type=int, default=32,
                        help="Samtools threads to read input BAM")

    parser.add_argument('--synthetic_proportion', type=float, default=0.25,
                        help="Target tumor proportion for synthetic pair data")

    parser.add_argument('--synthetic_coverage', type=int, default=None,
                        help="Target coverage for synthetic pair data")

    parser.add_argument('--contaminative_proportions', type=str, default=None,
                        help="contaminative_proportions, split by ','. ")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--cov_dir', type=str, default=None,
                        help="Directory of mosdepth coverage summary")

    # options for advanced users
    parser.add_argument('--tensor_sample_mode', type=str2bool, default=0,
                        help="EXPERIMENTAL: Add all tumor tensor and only sampling in tensor generation")

    parser.add_argument('--normal_coverage_proportion', type=float, default=1,
                        help="EXPERIMENTAL: Normal synthetic normal and tumor pair proportion")

    parser.add_argument('--dry_run', type=str2bool, default=0,
                        help="EXPERIMENTAL: Only print the synthetic log, debug only")


    # parser.add_argument('--normal_bam_fn', type=str, default=None,
    #                     help="Sorted normal BAM file input, required")
    #
    # parser.add_argument('--tumor_bam_fn', type=str, default=None,
    #                     help="Sorted tumor BAM file input, required")

    args = parser.parse_args()

    MixBin(args)


if __name__ == "__main__":
    main()
