import shlex
import os
import sys
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from shared.intervaltree.intervaltree import IntervalTree
from subprocess import run

from shared.utils import subprocess_popen

def select_hetero_snp_for_phasing(args):

    """
    Filter heterozygous snp variant for phasing, currently, we only filter snp variant with low quality socore as low
    quality variant contains more false positive variant that would lead to a larger minimum error correction loss.
    """
    qual_fn = args.qual_fn if args.qual_fn is not None else 'phase_qual'
    tumor_vcf_fn = args.tumor_vcf_fn
    normal_vcf_fn = args.normal_vcf_fn
    var_pct_full = args.var_pct_full
    contig_name = args.ctg_name
    output_folder = args.output_folder
    variant_dict = defaultdict(str)
    normal_qual_dict = defaultdict(int)
    tumor_qual_dict = defaultdict(int)
    found_qual_cut_off = False
    header = []

    #try to find the global quality cut off:
    f_qual = os.path.join(output_folder, qual_fn)
    if os.path.exists(f_qual):
        phase_qual_cut_off = float(open(f_qual, 'r').read().rstrip())
        found_qual_cut_off = True

    tumor_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (tumor_vcf_fn)))
    for row in tumor_unzip_process.stdout:
        row = row.rstrip()
        if row[0] == '#':
            header.append(row + '\n')
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        genotype = columns[9].split(':')[0].replace('|', '/')

        if len(ref_base) == 1 and len(alt_base) == 1:
            if genotype == '0/1' or genotype == '1/0':
                qual = float(columns[5])
                tumor_qual_dict[pos] = qual
                variant_dict[pos] = [ref_base, alt_base, qual, row]

    intersect_pos_set = defaultdict(int)
    intersect_pos_set = set()
    hetero_snp_not_found_in_tumor = 0
    hetero_snp_not_match_in_tumor = 0
    normal_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (normal_vcf_fn)))
    for row in normal_unzip_process.stdout:
        row = row.rstrip()
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        genotype = columns[9].split(':')[0].replace('|', '/')

        if len(ref_base) == 1 and len(alt_base) == 1:
            if genotype == '0/1' or genotype=='1/0':
                qual = float(columns[5])
                normal_qual_dict[pos] = qual
                if pos not in variant_dict:
                    hetero_snp_not_found_in_tumor += 1
                    continue
                tumor_ref_base , tumor_alt_base = variant_dict[pos][:2]
                if tumor_ref_base != ref_base or tumor_alt_base != alt_base:
                    hetero_snp_not_match_in_tumor += 1
                    continue
                intersect_pos_set.add(pos)
                # variant_dict[pos][-1] += ':' + ':'.join(row.split(':')[-2:])

    var_pct_full = 0.1
    # if found_qual_cut_off:
    #     remove_low_qual_list = [[k,v] for k,v in normal_qual_dict.items() if v < phase_qual_cut_off ]
    # else:
    normal_low_qual_set = set([item[0] for item in sorted(normal_qual_dict.items(), key=lambda x: x[1])[:int(var_pct_full * len(normal_qual_dict))]])
    tumor_low_qual_set = set([item[0] for item in sorted(tumor_qual_dict.items(), key=lambda x: x[1])[:int(var_pct_full * len(tumor_qual_dict))]])


    pass_variant_dict = defaultdict()
    low_qual_count = 0
    for pos in intersect_pos_set:
        if pos in normal_low_qual_set or pos in tumor_low_qual_set:
            low_qual_count += 1
            continue
        if pos in variant_dict:
            pass_variant_dict[pos] = variant_dict[pos][-1]

    print ('[INFO] Total heterozygous SNP positions selected: {}: {} not found:{} not match:{}, low_qual_count:{} normal:{} tumor:{}'.format(contig_name, len(pass_variant_dict), hetero_snp_not_found_in_tumor, hetero_snp_not_match_in_tumor, low_qual_count, len(normal_qual_dict), len(tumor_qual_dict)))

    if not os.path.exists(output_folder):
        return_code = run("mkdir -p {}".format(output_folder), shell=True)
    f = open(os.path.join(output_folder, '{}.vcf'.format(contig_name)), 'w')
    f.write(''.join(header))
    for key,row in sorted(pass_variant_dict.items(), key=lambda x: x[0]):
        f.write(row +'\n')
    f.close()


def main():
    parser = ArgumentParser(description="Select heterozygous snp candidates for WhatsHap phasing")

    parser.add_argument('--output_folder', type=str, default=None,
                        help="Path to directory that stores small bed region for raw alignment. (default: %(default)s)")

    parser.add_argument('--tumor_vcf_fn', type=str, default=None,
                        help="Path of the input vcf file. (default: %(default)s)")

    parser.add_argument('--normal_vcf_fn', type=str, default=None,
                        help="Path of the input vcf file. (default: %(default)s)")

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Default variant call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_pct_full', type=float, default=None,
                        help="Default reference call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--phase', action='store_false',
                        help="Only select hete candidates for phasing, default: True")

    parser.add_argument('--sampleName', type=str, default="",
                        help="Define the sample name to be shown in the VCF file, optional")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Input bam or sam have phasing info in HP tag, default: False")

    parser.add_argument('--split_bed_size', type=int, default=1000,
                        help="DEBUG: Default split bed size for parallel excution, default: %(default)s")

    parser.add_argument('--calling', type=int, default=0,
                        help="DEBUG: Path of the output folder, default: %(default)s")

    parser.add_argument('--realign_window_size', type=int, default=None,
                        help="DEBUG: The window size of read realignment, work with need_realignment")

    parser.add_argument('--split_region_size', type=int, default=40000000,
                        help="DEBUG: Vcf phasing split_region_size default: %(default)s")

    # options for internal process control
    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Output all alternative candidates path
    parser.add_argument('--all_alt_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Default chr prefix for contig name
    parser.add_argument('--chr_prefix', type=str, default='chr',
                        help=SUPPRESS)

    ## Default subsample depth for subsample bam file, 1000 means no subsampling
    parser.add_argument('--depth', type=int, default=1000,
                        help=SUPPRESS)

    ## Path of provided alternative file
    parser.add_argument('--alt_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Input the file that contains the quality cut-off for selecting low-quality pileup calls for phasing and full-alignment calling
    parser.add_argument('--qual_fn', type=str, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    select_hetero_snp_for_phasing(args)


if __name__ == "__main__":
    main()
