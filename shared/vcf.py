import shlex
import os

from textwrap import dedent
from subprocess import run
from collections import defaultdict

from shared.utils import subprocess_popen, Position as Position, file_path_from
import shared.param as param

caller_name = param.caller_name
version = param.version

vcf_header = dedent("""\
            ##fileformat=VCFv4.2
            ##source=ClairS
            ##{}_version={}
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Low quality variant">
            ##FILTER=<ID=RefCall,Description="Reference call">
            ##FILTER=<ID=Germline,Description="Germline variant">
            ##INFO=<ID=H,Number=0,Type=Flag,Description="Variant found only in one haplotype in the phased reads">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth in the tumor BAM">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the tumor BAM">
            ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth in the tumor BAM">
            ##FORMAT=<ID=NAF,Number=1,Type=Float,Description="Estimated allele frequency in the normal BAM">
            ##FORMAT=<ID=NDP,Number=1,Type=Integer,Description="Read depth in the normal BAM">
            ##FORMAT=<ID=NAD,Number=1,Type=Integer,Description="Alternate allele depth in the normal BAM">
            ##FORMAT=<ID=AU,Number=1,Type=Integer,Description="Count of A in the tumor BAM">
            ##FORMAT=<ID=CU,Number=1,Type=Integer,Description="Count of C in the tumor BAM">
            ##FORMAT=<ID=GU,Number=1,Type=Integer,Description="Count of G in the tumor BAM">
            ##FORMAT=<ID=TU,Number=1,Type=Integer,Description="Count of T in the tumor BAM">
            ##FORMAT=<ID=NAU,Number=1,Type=Integer,Description="Count of A in the normal BAM">
            ##FORMAT=<ID=NCU,Number=1,Type=Integer,Description="Count of C in the normal BAM">
            ##FORMAT=<ID=NGU,Number=1,Type=Integer,Description="Count of G in the normal BAM">
            ##FORMAT=<ID=NTU,Number=1,Type=Integer,Description="Count of T in the normal BAM">
            """.format(caller_name, version)
                    )


class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


class VcfWriter(object):
    def __init__(self,
                 vcf_fn,
                 ctg_name=None,
                 ref_fn=None,
                 sample_name="SAMPLE",
                 write_header=True,
                 header=None,
                 cmdline=None,
                 show_ref_calls=False):
        self.vcf_fn = vcf_fn
        self.show_ref_calls = show_ref_calls
        # make directory if not exist
        vcf_folder = os.path.dirname(self.vcf_fn)
        if not os.path.exists(vcf_folder):
            print("[INFO] Output VCF folder {} not found, create it".format(vcf_folder))
            return_code = run("mkdir -p {}".format(vcf_folder), shell=True)

        self.vcf_writer = open(self.vcf_fn, 'w')
        self.ref_fn = ref_fn
        self.ctg_name = ctg_name
        if ctg_name is not None:
            self.ctg_name_list = ctg_name.split(',') if ',' in ctg_name else [ctg_name]
        else:
            self.ctg_name_list = None
        self.sample_name = sample_name
        if write_header:
            self.write_header(ref_fn=ref_fn, header=header, cmdline=cmdline)

    def close(self):
        try:
            self.vcf_writer.close()
        except:
            pass

    def write_header(self, ctg_name=None, ref_fn=None, header=None, cmdline=None):
        header = vcf_header if header is None else header
        if cmdline is not None and cmdline != "":
            header_list = header.rstrip('\n').split('\n')
            insert_index = 3 if len(header_list) >= 3 else len(header_list) - 1
            header_list.insert(insert_index, "##cmdline={}".format(cmdline))
            header = "\n".join(header_list) + '\n'
        if self.ref_fn is not None:
            reference_index_file_path = file_path_from(self.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')
            with open(reference_index_file_path, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    if self.ctg_name_list is not None and contig_name not in self.ctg_name_list:
                        continue
                    header += "##contig=<ID=%s,length=%s>\n" % (contig_name, contig_size)

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (self.sample_name)

        self.vcf_writer.write(header)

    def write_row(self,
                  POS=None,
                  REF=None,
                  ALT=None,
                  QUAL=0,
                  GT='0/0',
                  DP=0,
                  AF=0,
                  AD=None,
                  NAD=None,
                  CHROM=None,
                  GQ=None,
                  ID='.',
                  FILTER=".",
                  INFO='.',
                  NAF=None,
                  TAF=None,
                  VT=None,
                  NDP=None,
                  TDP=None,
                  AU=None,
                  CU=None,
                  GU=None,
                  TU=None,
                  NAU=None,
                  NCU=None,
                  NGU=None,
                  NTU=None,
                  row_str=None):
        if row_str is not None:
            self.vcf_writer.write(row_str)
            return
        GQ = GQ if GQ else QUAL
        CHROM = CHROM if CHROM else self.ctg_name
        if not self.show_ref_calls and (GT == "0/0" or GT == "./."):
            return
        FORMAT = "GT:GQ:DP:AF"
        FORMAT_V = "%s:%.4f:%d:%.4f" % (GT, GQ, DP, AF)
        basic_vcf_format = "%s\t%d\t%s\t%s\t%s\t%.4f\t%s\t%s" % (
            CHROM,
            int(POS),
            ID,
            REF,
            ALT,
            QUAL,
            FILTER,
            INFO
        )
        if AD is not None and AD != "":
            FORMAT += ":AD"
            FORMAT_V += ":%d" % (AD)
        if NAF is not None:
            FORMAT += ":NAF"
            FORMAT_V += ":%.4f" % (NAF)
        if TAF is not None:
            FORMAT += ":TAF"
            FORMAT_V += ":%.4f" % (TAF)
        if NDP is not None:
            FORMAT += ":NDP"
            FORMAT_V += ":%d" % (NDP)
        if NAD is not None:
            FORMAT += ":NAD"
            FORMAT_V += ":%d" % (NAD)
        if TDP is not None:
            FORMAT += ":TDP"
            FORMAT_V += ":%d" % (TDP)
        if AU is not None and CU is not None and GU is not None and TU is not None:
            FORMAT += ":AU:CU:GU:TU"
            FORMAT_V += ":%d:%d:%d:%d" % (AU, CU, GU, TU)
        if NAU is not None and NCU is not None and NGU is not None and NTU is not None:
            FORMAT += ":NAU:NCU:NGU:NTU"
            FORMAT_V += ":%d:%d:%d:%d" % (NAU, NCU, NGU, NTU)

        if VT is not None:
            FORMAT += ":VT"
            FORMAT_V += ":%s" % (VT)
        vcf_format = '\t'.join([basic_vcf_format, FORMAT, FORMAT_V]) + "\n"

        self.vcf_writer.write(vcf_format)


class VcfReader(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 ctg_start=None,
                 ctg_end=None,
                 is_var_format=False,
                 is_happy_format=False,
                 is_fp=None,
                 show_ref=True,
                 direct_open=False,
                 keep_row_str=False,
                 skip_genotype=False,
                 filter_tag=None,
                 naf_filter=None,
                 taf_filter=None,
                 save_header=False,
                 min_qual=None,
                 max_qual=None,
                 discard_indel=False,
                 keep_af=False):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.ctg_start = ctg_start
        self.ctg_end = ctg_end
        self.variant_dict = defaultdict(Position)
        self.is_var_format = is_var_format
        self.is_happy_format = is_happy_format
        self.is_fp = is_fp
        self.show_ref = show_ref
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.skip_genotype = skip_genotype
        self.filter_tag = filter_tag  # PASS;HighConf PASS;MedConf in hcc1395
        self.naf_filter = naf_filter
        self.taf_filter = taf_filter
        self.header = ""
        self.save_header = save_header
        self.discard_indel = discard_indel
        self.min_qual = min_qual
        self.max_qual = max_qual
        self.keep_af = keep_af

    def read_vcf(self):
        is_ctg_region_provided = self.ctg_start is not None and self.ctg_end is not None

        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        header_last_column = []
        if self.direct_open:
            vcf_fp = open(self.vcf_fn)
            vcf_fo = vcf_fp
        else:
            vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
            vcf_fo = vcf_fp.stdout
        for row in vcf_fo:
            columns = row.strip().split()
            if columns[0][0] == "#":
                if self.save_header:
                    self.header += row
                header_last_column = columns
                continue

            tumor_in_last = True if len(header_last_column) and header_last_column[
                -1].rstrip().lower() == "tumor" else False
            # position in vcf is 1-based
            chromosome, position = columns[0], columns[1]
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue
            if is_ctg_region_provided and not (self.ctg_start <= int(position) <= self.ctg_end):
                continue

            # filter NAF
            if self.naf_filter is not None:
                naf_index = columns[8].split(':').index("NAF")
                if float(columns[9].split(':')[naf_index]) >= self.naf_filter:
                    continue
            if self.taf_filter is not None:
                taf_index = columns[8].split(':').index("AF")
                naf_index = columns[8].split(':').index("NAF")
                if float(columns[9].split(':')[taf_index]) * self.taf_filter < float(columns[9].split(':')[naf_index]):
                    continue

            FILTER = columns[6] if len(columns) >= 7 else None
            if self.filter_tag is not None:
                filter_list = self.filter_tag.split(',')
                if sum([1 if filter == FILTER else 0 for filter in filter_list]) == 0:
                    continue
            self.is_var_format = True if columns[2][0] in 'ACGT' else False
            self.is_var_format = False
            if self.is_var_format:
                reference, alternate = columns[2], columns[3]
                genotype_1 = int(columns[4])
                genotype_2 = int(columns[5])
            else:
                reference, alternate, last_column = columns[3], columns[4], columns[-1]

                if self.discard_indel and (len(reference) > 1 or len(alternate) > 1):
                    continue

                try:
                    qual = columns[5] if len(columns) > 5 else None

                    if self.min_qual is not None and float(qual) < self.min_qual:
                        continue

                    if self.max_qual is not None and float(qual) > self.max_qual:
                        continue
                except:
                    qual = None

                last_column = last_column if not tumor_in_last else columns[-2]
                if self.is_happy_format and self.is_fp:
                    last_column = columns[10]
                if self.is_happy_format and not self.is_fp:
                    last_column = columns[9]
                genotype = last_column.split(":")[0].replace("/", "|").replace(".", "0").split("|")
                try:
                    genotype_1, genotype_2 = genotype

                    if int(genotype_1) > int(genotype_2):
                        genotype_1, genotype_2 = genotype_2, genotype_1

                    # remove * to guarentee vcf match
                    if '*' in alternate:
                        alternate = alternate.split(',')
                        if int(genotype_1) + int(genotype_2) != 3 or len(alternate) != 2:
                            print('error with variant representation')
                            continue
                        alternate = ''.join([alt_base for alt_base in alternate if alt_base != '*'])
                        # * always have a genotype 1/2

                        genotype_1, genotype_2 = '0', '1'
                except:
                    genotype_1 = -1
                    genotype_2 = -1
            if self.keep_af:
                tag_list = columns[8].split(':')
                if 'AF' in tag_list or 'VAF' in tag_list:
                    taf_index = tag_list.index('AF') if 'AF' in tag_list else tag_list.index('VAF')
                    taf = float(columns[9].split(':')[taf_index])
                else:
                    taf = None
            else:
                taf = None
            position = int(position)
            have_extra_infos = 'VT' in row

            if genotype_1 == "0" and genotype_2 == "0" and not self.show_ref and not self.skip_genotype:
                continue
            extra_infos = columns[-1].split(':')[-1] if have_extra_infos else ''
            row_str = row if self.keep_row_str else False
            key = (chromosome, position) if self.ctg_name is None else position

            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              ref_base=reference,
                                              alt_base=alternate,
                                              genotype1=int(genotype_1),
                                              genotype2=int(genotype_2),
                                              qual=qual,
                                              row_str=row_str,
                                              af=taf,
                                              filter=FILTER,
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
