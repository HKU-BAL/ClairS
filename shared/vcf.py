import sys
import shlex

from subprocess import PIPE
from argparse import ArgumentParser
from collections import defaultdict
from shared.utils import subprocess_popen, Position as Position

class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()

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

        if self.vcf_fn is None:
            return

        header_last_column = []
        vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
        for row in vcf_fp.stdout:
            columns = row.strip().split()
            if columns[0][0] == "#":
                header_last_column = columns
                continue

            tumor_first = header_last_column[-1].rstrip().lower() != "tumor" if len(header_last_column) else True
            # position in vcf is 1-based
            chromosome, position = columns[0], columns[1]
            if chromosome != self.ctg_name:
                continue
            if is_ctg_region_provided and not (self.ctg_start <= int(position) <= self.ctg_end):
                continue
            if self.is_var_format:
                reference, alternate = columns[2], columns[3]
                genotype_1 = int(columns[4])
                genotype_2 = int(columns[5])
            else:
                reference, alternate, last_column = columns[3], columns[4], columns[-1]
            # normal GetTruth
                last_column = last_column if not tumor_first else columns[-2]
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
            self.variant_dict[position] = Position(pos=position,
                                                    ref_base=reference,
                                                   alt_base=alternate,
                                                   genotype1=int(genotype_1),
                                                   genotype2=int(genotype_2))
    def get_alt_info(self, pos, extra_info=""):
        pos = int(pos)
        if pos not in self.variant_dict:
            return ""
        ref_base = self.variant_dict[pos].reference_bases
        alt_base = ','.join(self.variant_dict[pos].alternate_bases)
        gentoype_str = '/'.join([str(g) for g in self.variant_dict[pos].genotype])
        return extra_info + '_' + ref_base + '_' + alt_base + '_' + gentoype_str