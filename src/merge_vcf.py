import os
import subprocess
import shlex
from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict
import sys
sys.path.insert(0, "/autofs/bal36/zxzheng/somatic/Clair-somatic")
import shlex
from shared.vcf import VcfReader

from shared.utils import log_error, log_warning, file_path_from, str2bool
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def output_header(output_fn, reference_file_path, sample_name='SAMPLE'):
    output_file = open(output_fn, "w")
    from textwrap import dedent
    output_file.write(dedent("""\
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FILTER=<ID=LowQual,Description="Low quality variant">
        ##FILTER=<ID=RefCall,Description="Reference call">
        ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
        ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">"""
                  ) + '\n')

    if reference_file_path is not None:
        reference_index_file_path = file_path_from(reference_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
        with open(reference_index_file_path, "r") as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")
                contig_name, contig_size = columns[0], columns[1]
                output_file.write(("##contig=<ID=%s,length=%s>" % (contig_name, contig_size) + '\n'))

    output_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name))
    output_file.close()

def print_calling_step(output_fn=""):

    merge_output = os.path.join(os.path.dirname(output_fn), 'merge_output.vcf.gz')
    pileup_output = os.path.join(os.path.dirname(output_fn), 'pileup.vcf.gz')

    # print (log_warning("[WARNING] Copying pileup.vcf.gz to {}".format(merge_output)))
    # subprocess.run('cp {} {}'.format(pileup_output, merge_output), shell=True, stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE)

def check_header_in_gvcf(header, contigs_list):
    # Only output the contigs processed to be consistent with GATK
    # Contig format: ##contig=<ID=%s,length=%s>

    update_header = []
    for row_id, row in enumerate(header):
        if row.startswith("##contig="):
            contig = row.split(',')[0].split('=')[2]
            if contig not in contigs_list:
                continue
        update_header.append(row)

    return update_header

def merge_vcf(args):

    compress_vcf = args.compress_vcf

    input_vcf_reader = VcfReader(vcf_fn=args.full_alignment_vcf_fn, ctg_name=None, show_ref=False, keep_row_str=True,
                                 skip_genotype=True,
                                 filter_tag="PASS",
                                 keep_af=True)  # , naf_filter=0.03, taf_filter=0.25)
    input_vcf_reader.read_vcf()
    fa_input_variant_dict = input_vcf_reader.variant_dict

    # if args.pileup_vcf_fn is not None:
    #     from shared.vcf import VcfReader
    #     p_input_vcf_reader = VcfReader(vcf_fn=args.pileup_vcf_fn, ctg_name=None, show_ref=False, keep_row_str=True,
    #                                  skip_genotype=True,
    #                                  filter_tag="PASS",
    #                                  keep_af=True)  # , naf_filter=0.03, taf_filter=0.25)
    #     p_input_vcf_reader.read_vcf()
    #     p_input_variant_dict = p_input_vcf_reader.variant_dict
    #
    # else:

    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    filter_count = 0
    af_filter_count = 0
    f = open(args.pileup_vcf_fn)
    for row in f:

        if row[0] == '#':
            if row not in header:
                header.append(row)
            continue
        row_count += 1
        # use the first vcf header
        columns = row.strip().split()
        ctg_name, pos = columns[0], columns[1]
        # print(columns)
        qual = float(columns[5])
        filter = columns[6]
        if filter != 'PASS':
            continue
        if args.qual is not None and qual <= args.qual:

            if (ctg_name, int(pos)) not in fa_input_variant_dict:
                filter_count += 1
                continue


        if args.af is not None:
            tag_list = columns[8].split(':')
            taf_index = tag_list.index('AF') if 'AF' in tag_list else tag_list.index('VAF')
            af = float(columns[9].split(':')[taf_index])
            if af <= args.af and (ctg_name, int(pos)) not in fa_input_variant_dict:
                af_filter_count += 1
                continue

        if args.qual is not None and (ctg_name, int(pos)) in fa_input_variant_dict:
            columns[5] = str((qual + float(fa_input_variant_dict[(ctg_name, int(pos))].qual)) / 2)
            # columns[5] = str(input_variant_dict[(ctg_name, int(pos))].qual)
            row = '\t'.join(columns) + '\n'

        contig_dict[ctg_name][int(pos)] = row
        no_vcf_output = False

    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, please check the setting"))
    if no_vcf_output:
        print(log_warning("[WARNING] No variant found, please check the setting"))

    print("Filter count", filter_count)
    if args.af is not None:
        print("AF filter count", af_filter_count)
    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(args.output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])
    f.close()

    if compress_vcf:
        compress_index_vcf(args.output_fn)

def sort_bed_from_stdin(args):
    """
    Sort vcf file according to variants start position and contig name.
    """

    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    for row in stdin:
        row_count += 1
        if row[0] == '#':
            if row not in header:
                header.append(row)
            continue
        # use the first vcf header
        columns = row.strip().split()
        ctg_name, bed_start, bed_end = columns[:3]
        contig_dict[ctg_name][(int(bed_start), int(bed_end))] = row
        no_bed_output = False
    if row_count == 0:
        print(log_warning("[WARNING] No BED file found, please check the setting"))
    if no_bed_output:
        print(log_warning("[WARNING] No BED found, please check the setting"))

    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(args.output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])


def sort_vcf_from(args):
    """
    Sort vcf file from providing vcf filename prefix.
    """
    output_fn = args.output_fn
    input_dir = args.input_dir
    vcf_fn_prefix = args.vcf_fn_prefix
    vcf_fn_suffix = args.vcf_fn_suffix
    sample_name = args.sampleName
    ref_fn = args.ref_fn
    contigs_fn = args.contigs_fn
    compress_vcf = args.compress_vcf

def mergeNonVariant(args):
    '''
    merge the variant calls and non-variants

    '''
    gvcf_generator = gvcfGenerator(ref_path=args.ref_fn, samtools=args.samtools)
    raw_gvcf_path = args.non_var_gvcf_fn
    raw_vcf_path = args.output_fn
    
    if (args.gvcf_fn == None):
        save_path = args.call_fn.split('.')[0] + '.g.vcf'
    else:
        save_path = args.gvcf_fn
        logging.info("[INFO] Merge variants and non-variants to GVCF")
        gvcf_generator.mergeCalls(raw_vcf_path, raw_gvcf_path, save_path, args.sampleName, args.ctgName, args.ctgStart,
                                  args.ctgEnd)
    pass


def main():
    parser = ArgumentParser(description="Generate 1-based variant candidates using alignments")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Path to the pileup vcf file")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Path to the full alignment vcf file")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--non_var_gvcf_fn', type=str, default=None,
                        help='Path to the non-variant GVCF')

    parser.add_argument('--gvcf_fn', type=str, default=None,
                        help='Filename of the GVCF output')

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Filename of the merged output")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--bed_fn_prefix', type=str, default=None,
                        help="Process variant only in the provided regions prefix")

    parser.add_argument('--qual', type=int, default=2,
                        help="If set, variants with >$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file")

    parser.add_argument('--samtools', type=str, default='samtools',
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--bed_format', action='store_true',
                        help="Only work for gvcf file, reduce hard disk space")

    parser.add_argument('--filter_vcf_fn', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--qual', type=float, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--af', type=float, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--compress_vcf', type=str2bool, default=True,
                        help="Compress and index output VCF")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    
    # realignment region merge
    if args.platform == 'ilmn':
        MergeVcf_illumina(args=args)
    else:
        MergeVcf(args=args)
    
    if (args.gvcf):
        mergeNonVariant(args)


if __name__ == "__main__":
    main()
