import sys
import shlex
import os

from textwrap import dedent
from subprocess import PIPE
from argparse import ArgumentParser
from collections import defaultdict
from shared.utils import subprocess_popen, Position as Position, file_path_from

class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


class VcfWriter(object):
    def __init__(self, vcf_fn, ctg_name=None, ref_fn = None, write_header=True, show_ref_calls=False):
        self.vcf_fn = vcf_fn
        self.show_ref_calls = show_ref_calls
        self.vcf_writer = open(self.vcf_fn, 'w')
        self.ref_fn = ref_fn
        self.ctg_name = ctg_name
        if write_header:
            self.write_header(ref_fn=ref_fn)

    def close(self):
        try:
            self.vcf_writer.close()
        except:
            pass

    def write_header(self, ctg_name=None, ref_fn=None, sample_name="SAMPLE"):
        header = dedent("""\
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
                    ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">
                    """
               )

        if self.ref_fn is not None:
            reference_index_file_path = file_path_from(self.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')
            with open(reference_index_file_path, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    if self.ctg_name is not None and contig_name != self.ctg_name:
                        continue
                    header += "##contig=<ID=%s,length=%s>\n" % (contig_name, contig_size)

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (sample_name)

        self.vcf_writer.write(header)

    def write_row(self, POS, REF, ALT, QUAL, GT, DP, AF, CHROM=None, GQ=None, ID='.', FILTER=".", INFO='.', NAF=None, TAF=None, VT=None):
        GQ = GQ if GQ else QUAL
        CHROM = CHROM if CHROM else self.ctg_name
        if not self.show_ref_calls and (GT == "0/0" or GT == "./."):
            return
        FORMAT = "GT:GQ:DP:AF"
        FORMAT_V = "%s:%d:%d:%.4f" % (GT, GQ, DP, AF)
        basic_vcf_format = "%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s" % (
            CHROM,
            int(POS),
            ID,
            REF,
            ALT,
            QUAL,
            FILTER,
            INFO
            )
        if NAF is not None:
            FORMAT += ":NAF"
            FORMAT_V += ":%.4f" % (NAF)
        if TAF is not None:
            FORMAT += ":TAF"
            FORMAT_V += ":%.4f" % (TAF)
        if VT is not None:
            FORMAT += ":VT"
            FORMAT_V += ":%s" % (VT)

        vcf_format = '\t'.join([basic_vcf_format, FORMAT, FORMAT_V]) + "\n"

        self.vcf_writer.write(vcf_format)

class VcfReader(object):
    def __init__(self, vcf_fn, ctg_name, ctg_start=None, ctg_end=None, is_var_format=True, is_happy_format=False, is_fp=None):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.ctg_start = ctg_start
        self.ctg_end = ctg_end
        self.variant_dict = defaultdict(Position)
        self.is_var_format = is_var_format
        self.is_happy_format = is_happy_format
        self.is_fp = is_fp
    def read_vcf(self):
        is_ctg_region_provided = self.ctg_start is not None and self.ctg_end is not None

        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        header_last_column = []
        vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
        for row in vcf_fp.stdout:
            columns = row.strip().split()
            if columns[0][0] == "#":
                header_last_column = columns
                continue

            tumor_in_last = True if len(header_last_column) and header_last_column[-1].rstrip().lower() == "tumor" else False
            # position in vcf is 1-based
            chromosome, position = columns[0], columns[1]
            if chromosome != self.ctg_name:
                continue
            if is_ctg_region_provided and not (self.ctg_start <= int(position) <= self.ctg_end):
                continue
            self.is_var_format = True if columns[2][0] in 'ACGT' else False
            if self.is_var_format:
                reference, alternate = columns[2], columns[3]
                genotype_1 = int(columns[4])
                genotype_2 = int(columns[5])
            else:
                reference, alternate, last_column = columns[3], columns[4], columns[-1]
            # normal GetTruth
                last_column = last_column if not tumor_in_last else columns[-2]
                if self.is_happy_format and self.is_fp:
                    last_column = columns[10]
                if self.is_happy_format and not self.is_fp:
                    last_column = columns[9]
                genotype = last_column.split(":")[0].replace("/", "|").replace(".", "0").split("|")
                try:
                    genotype_1, genotype_2 = genotype

                    # 1000 Genome GetTruth (format problem) (no genotype is given)
                    if int(genotype_1) > int(genotype_2):
                        genotype_1, genotype_2 = genotype_2, genotype_1

                    #remove * to guarentee vcf match
                    if '*' in alternate:
                        alternate = alternate.split(',')
                        if int(genotype_1) + int(genotype_2) != 3 or len(alternate) != 2:
                            print ('error with variant representation')
                            continue
                        alternate = ''.join([alt_base for alt_base in alternate if alt_base != '*'])
                        # * always have a genotype 1/2

                        genotype_1, genotype_2 = '0', '1'
                except:
                    genotype_1 = -1
                    genotype_2 = -1
            position = int(position)
            have_extra_infos = 'VT' in row

            extra_infos = columns[-1].split(':')[-1] if have_extra_infos else ''
            self.variant_dict[position] = Position(pos=position,
                                                    ref_base=reference,
                                                   alt_base=alternate,
                                                   genotype1=int(genotype_1),
                                                   genotype2=int(genotype_2),
                                                   extra_infos=extra_infos)
    def get_alt_info(self, pos, extra_info=""):
        pos = int(pos)
        if pos not in self.variant_dict:
            return ""
        ref_base = self.variant_dict[pos].reference_bases
        alt_base = ','.join(self.variant_dict[pos].alternate_bases)
        gentoype_str = '/'.join([str(g) for g in self.variant_dict[pos].genotype])
        extra_info = self.variant_dict[pos].extra_infos if self.variant_dict[pos].extra_infos != "" else extra_info
        return extra_info + '_' + ref_base + '_' + alt_base + '_' + gentoype_str