import os
import subprocess
import shlex
from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict
import sys
import os
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

def update_haploid_sensitive_genotype(columns):
    INFO = columns[9].split(':')
    genotype_string = INFO[0].replace('|', '/')
    ref_base, alt_base = columns[3], columns[4]
    is_multi = ',' in alt_base

    if is_multi:
        return ""

    if genotype_string in ('0/1','1/0','1/1'):
        genotype = ['1']
    else:
        genotype = ['0']
    # update genotype
    columns[9] = ':'.join(genotype + INFO[1:])
    row = '\t'.join(columns) + '\n'
    return row

def MarkLowQual(row, quality_score_for_pass, qual):
    if row == '':
        return row

    if quality_score_for_pass and qual <= quality_score_for_pass:
        row = row.split("\t")
        row[6] = "LowQual"
        return '\t'.join(row)
    return row


def MergeVcf(args):
    """
    Merge pileup and full alignment vcf output. We merge the low quality score pileup candidates
    recalled by full-alignment model with high quality score pileup output.
    """

    output_fn = args.output_fn
    full_alignment_vcf_fn = args.full_alignment_vcf_fn
    pileup_vcf_fn = args.pileup_vcf_fn  # true vcf var
    contig_name = args.ctgName
    QUAL = args.qual
    is_haploid_precise_mode_enabled = args.haploid_precise
    is_haploid_sensitive_mode_enabled = args.haploid_sensitive
    print_ref = args.print_ref_calls
    full_alignment_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (full_alignment_vcf_fn)))

    full_alignment_output = []
    full_alignment_output_set = set()
    header = []

    for row in full_alignment_vcf_unzip_process.stdout:
        if row[0] == '#':
            header.append(row)
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue
        pos = int(columns[1])
        qual = float(columns[5])
        ref_base, alt_base = columns[3], columns[4]
        is_reference = (alt_base == "." or ref_base == alt_base)

        full_alignment_output_set.add((ctg_name, pos))

        if is_haploid_precise_mode_enabled:
            row = update_haploid_precise_genotype(columns)
        if is_haploid_sensitive_mode_enabled:
            row = update_haploid_sensitive_genotype(columns)

        if not is_reference:
                row = MarkLowQual(row, QUAL, qual)
                full_alignment_output.append((pos, row))

        elif print_ref:
            full_alignment_output.append((pos, row))

    full_alignment_vcf_unzip_process.stdout.close()
    full_alignment_vcf_unzip_process.wait()

    pileup_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (pileup_vcf_fn)))

    output_file = open(output_fn, 'w')
    output_file.write(''.join(header))
    
    def pileup_vcf_generator_from(pileup_vcf_unzip_process):
        pileup_row_count = 0
        for row in pileup_vcf_unzip_process.stdout:
            if row[0] == '#':
                continue

            columns = row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            qual = float(columns[5])
            ref_base, alt_base = columns[3], columns[4]
            is_reference = (alt_base == "." or ref_base == alt_base)

            if (ctg_name, pos) in full_alignment_output_set:
                continue

            if is_haploid_precise_mode_enabled:
                row = update_haploid_precise_genotype(columns)
            if is_haploid_sensitive_mode_enabled:
                row = update_haploid_sensitive_genotype(columns)

            if not is_reference:
                row = MarkLowQual(row, QUAL, qual)
                pileup_row_count += 1
                yield (pos, row)
            elif print_ref:
                pileup_row_count += 1
                yield (pos, row)

        logging.info('[INFO] Pileup variants processed in {}: {}'.format(contig_name, pileup_row_count))

    pileup_vcf_generator = pileup_vcf_generator_from(pileup_vcf_unzip_process=pileup_vcf_unzip_process)
    full_alignment_vcf_generator = iter(full_alignment_output)
    for vcf_infos in heapq.merge(full_alignment_vcf_generator, pileup_vcf_generator):
        if len(vcf_infos) != 2:
            continue
        pos, row = vcf_infos
        output_file.write(row)
    
    logging.info('[INFO] Full-alignment variants processed in {}: {}'.format(contig_name, len(full_alignment_output)))

    pileup_vcf_unzip_process.stdout.close()
    pileup_vcf_unzip_process.wait()
    output_file.close()

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

    parser.add_argument('--print_ref_calls', type=str2bool, default=False,
                        help="Show reference calls (0/0) in vcf file output")

    # options for advanced users
    parser.add_argument('--haploid_precise', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant")

    parser.add_argument('--haploid_sensitive', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant")

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
