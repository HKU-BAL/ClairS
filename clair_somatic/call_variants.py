import sys
import os
import math
# import tables
# import tensorflow as tf
# import numpy as np
import logging
import shlex
from time import time
from argparse import ArgumentParser, SUPPRESS
from threading import Thread
from subprocess import run
from math import log, e
from collections import namedtuple
from shared.vcf import VcfWriter

from clair_somatic.task.gt21 import (
    GT21_Type, gt21_enum_from_label,
    HOMO_SNP_GT21, HOMO_SNP_LABELS,
    HETERO_SNP_GT21, HETERO_SNP_LABELS, GT21_LABELS, partial_label_from, mix_two_partial_labels
)
# import clair_somatic.utils as utils
from clair_somatic.task.genotype import Genotype, genotype_string_from, genotype_enum_from, genotype_enum_for_task
from shared.utils import IUPAC_base_to_ACGT_base_dict as BASE2ACGT, BASIC_BASES, str2bool, file_path_from, log_error, log_warning, subprocess_popen
from clair_somatic.task.variant_length import VariantLength
import shared.param as param
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.basicConfig(format='%(message)s', level=logging.INFO)
minimum_variant_length_that_need_infer = VariantLength.max
maximum_variant_length_that_need_infer = 50
ACGT = 'ACGT'
Phred_Trans = (-10 * log(e, 10))

OutputConfig = namedtuple('OutputConfig', [
    'is_show_reference',
    'is_show_germline',
    'is_haploid_precise_mode_enabled',
    'is_haploid_sensitive_mode_enabled',
    'is_output_for_ensemble',
    'quality_score_for_pass',
    'tensor_fn',
    'input_probabilities',
    'add_indel_length',
    'gvcf',
    'pileup'
])
OutputUtilities = namedtuple('OutputUtilities', [
    'print_debug_message',
    'output',
    'output_header',
    'close_opened_files',
    'gen_output_file'
])


def homo_SNP_bases_from(homo_SNP_probabilities):
    output_bases = HOMO_SNP_LABELS[np.argmax(homo_SNP_probabilities)]
    return output_bases[0], output_bases[1]


def hetero_SNP_bases_from(hetero_SNP_probabilities):
    output_bases = HETERO_SNP_LABELS[np.argmax(hetero_SNP_probabilities)]
    return output_bases[0], output_bases[1]


def filtration_value_from(quality_score_for_pass, quality_score, is_reference=False, is_germline=False):
    """
    filter qual if set quality score, variant quliaty lower than specific quality socre will be
    marked ass LowQual otherwise PASS. Default there is no quality score cut off.

    """
    if is_reference:
        return 'RefCall'
    elif is_germline:
        return 'Germline'
    if quality_score_for_pass is None:
        return "PASS"
    if quality_score >= quality_score_for_pass:
        return "PASS"

    return "LowQual"


def insertion_bases_using_alt_info_from(
        alt_info_dict,
        propose_insertion_length=None,
        minimum_insertion_length=1,
        maximum_insertion_length=maximum_variant_length_that_need_infer,
        insertion_bases_to_ignore="",
        return_multi=False
):
    """
    get insertion base using altnertive information in bam alignment file.
    alt_info_dict: dictionary (XID: count), include snp, insertion, deletion type and read count.
    propose_insertion_length: if set, only return insertion length which match propose insertion length.
    minimum_insertion_length: if set, only return insertion length which is larger than specific insertion length.
    maximum_insertion_length: if set, only return insertion which is shorter than specific insertion length,
    we will always only return insertion length shorter than 50bp by default.
    insertion_bases_to_ignore: for multi alleic insertion variants, set the insertion bases to be ignored.
    """

    if propose_insertion_length:
        propose_insertion_length += 1  # include reference base
    if not len(alt_info_dict): return ""
    insertion_bases_dict = {}
    propose_insertion_bases_dict = {}
    for raw_key, items in alt_info_dict.items():
        if raw_key[0] != 'I': continue
        key = raw_key[1:]  # remove first cigar +-X and reference_base
        if propose_insertion_length and len(key) == propose_insertion_length and key != insertion_bases_to_ignore:
            propose_insertion_bases_dict[key] = items
        elif minimum_insertion_length <= len(key) <= maximum_insertion_length and key != insertion_bases_to_ignore:
            insertion_bases_dict[key] = items

    if propose_insertion_length and len(propose_insertion_bases_dict):
        return max(propose_insertion_bases_dict, key=propose_insertion_bases_dict.get) if len(
            propose_insertion_bases_dict) > 0 else ""
    if return_multi:
        insertion_bases_list = list(insertion_bases_dict.items())
        insertion_bases_list = [item[0] for item in sorted(insertion_bases_list, key=lambda x: x[1])[::-1]]
        return insertion_bases_list[:2] if len(insertion_bases_list) else ""

    return max(insertion_bases_dict, key=insertion_bases_dict.get) if len(insertion_bases_dict) > 0 else ""


def deletion_bases_using_alt_info_from(
        alt_info_dict,
        propose_deletion_length=None,
        minimum_deletion_length=1,
        maximum_deletion_length=maximum_variant_length_that_need_infer,
        deletion_bases_to_ignore="",
        return_multi=False,

):
    """
    get deletion base using altnertive information in bam alignment file.
    alt_info_dict: dictionary (XID: count), include snp, insertion, deletion type and read count.
    propose_deletion_length: if set, only return deletion length which match propose deletion length.
    minimum_deletion_length: if set, only return deletion length which is larger than specific deletion length.
    maximum_deletion_length: if set, only return deletion which is shorter than specific deletion length,
    we will always only return deletion length shorter than 50bp by default.
    deletion_bases_to_ignore: for multi alleic variants, set the deletion bases to be ignored.
    """

    if not len(alt_info_dict): return ""

    deletion_bases_dict = {}
    propose_deletion_bases_dict = {}
    for raw_key, items in alt_info_dict.items():
        if raw_key[0] != 'D': continue
        key = raw_key[1:]  # remove first cigar +-X
        if propose_deletion_length and len(key) == propose_deletion_length and key != deletion_bases_to_ignore:
            propose_deletion_bases_dict[key] = items

        elif minimum_deletion_length <= len(key) <= maximum_deletion_length and key != deletion_bases_to_ignore:
            deletion_bases_dict[key] = items

    if propose_deletion_length and len(propose_deletion_bases_dict):
        return max(propose_deletion_bases_dict, key=propose_deletion_bases_dict.get) if len(
            propose_deletion_bases_dict) > 0 else ""

    if return_multi:
        deletion_bases_list = list(deletion_bases_dict.items())
        deletion_bases_list = [item[0] for item in sorted(deletion_bases_list, key=lambda x: x[1])[::-1]]
        if len(deletion_bases_list) <= 1: return ""
        return [deletion_bases_list[0], deletion_bases_list[1]] if len(deletion_bases_list[0]) > len(
            deletion_bases_list[1]) else [deletion_bases_list[1], deletion_bases_list[0]]
    return max(deletion_bases_dict, key=deletion_bases_dict.get) if len(deletion_bases_dict) > 0 else ""


# def call_variants(args):
#     call_variants_with_probabilities_input(args=args)

def output_utilties_from(
        sample_name,
        is_output_for_ensemble,
        reference_file_path,
        output_file_path,
        output_probabilities
):
    def gen_output_file():
        global output_file
        if not output_probabilities:
            # output_file = VcfWriter(vcf_fn=output_file_path, ref_fn=reference_file_path, )
            output_file = open(output_file_path, "w")

    def output(string_value):
        global output_file
        print(string_value, file=output_file)

    def print_debug_message(
            chromosome,
            position,
            gt21_probabilities,
            genotype_probabilities,
            variant_length_probabilities_1,
            variant_length_probabilities_2,
            extra_infomation_string=""
    ):
        output("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            chromosome,
            position,
            ["{:0.8f}".format(x) for x in gt21_probabilities],
            ["{:0.8f}".format(x) for x in genotype_probabilities],
            ["{:0.8f}".format(x) for x in variant_length_probabilities_1],
            ["{:0.8f}".format(x) for x in variant_length_probabilities_2],
            extra_infomation_string
        ))

    def close_opened_files():
        output_file.close()

    def output_header():
        if is_output_for_ensemble:
            return

        from textwrap import dedent
        output(dedent("""\
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
                      ))

        if reference_file_path is not None:
            reference_index_file_path = file_path_from(reference_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
            with open(reference_index_file_path, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    output("##contig=<ID=%s,length=%s>" % (contig_name, contig_size))

        output('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name))

    return OutputUtilities(
        print_debug_message,
        output,
        output_header,
        close_opened_files,
        gen_output_file
    )


def homo_Ins_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2, extra_probability):
    return [(
        i,
        variant_length_probabilities_1[i + VariantLength.index_offset] *
        variant_length_probabilities_2[i + VariantLength.index_offset] * extra_probability
    ) for i in range(1, VariantLength.max + 1)]


def hetero_Ins_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
    return [(
        i,
        variant_length_probabilities_1[0 + VariantLength.index_offset] *
        variant_length_probabilities_2[i + VariantLength.index_offset],

    ) for i in range(1, VariantLength.max + 1)]


def hetero_InsIns_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2, extra_probability):
    probabilities = []
    for i in range(1, VariantLength.max + 1):
        for j in range(i, VariantLength.max + 1):
            # note: one kind of InsIns is same # of insertion bases but different kind of ACGT
            probabilities.append((
                (i, j),
                variant_length_probabilities_1[i + VariantLength.index_offset] *
                variant_length_probabilities_2[j + VariantLength.index_offset] * extra_probability
            ))
    return probabilities


def homo_Del_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2, extra_probability):
    return [(
        i,
        variant_length_probabilities_1[-i + VariantLength.index_offset] *
        variant_length_probabilities_2[-i + VariantLength.index_offset] * extra_probability
    ) for i in range(1, VariantLength.max + 1)]


def hetero_Del_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
    return [(
        i,
        variant_length_probabilities_1[-i + VariantLength.index_offset] *
        variant_length_probabilities_2[0 + VariantLength.index_offset],
    ) for i in range(1, VariantLength.max + 1)]


def hetero_DelDel_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2, extra_probability):
    probabilities = []
    for i in range(1, VariantLength.max + 1):
        for j in range(1, VariantLength.max + 1):
            if i == j and i != VariantLength.index_offset and j != VariantLength.index_offset:
                continue
            probabilities.append((
                (i, j) if i < j else (j, i),
                variant_length_probabilities_1[-i + VariantLength.index_offset] *
                variant_length_probabilities_2[-j + VariantLength.index_offset] * extra_probability
            ))
    return probabilities


def hetero_InsDel_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2, extra_probability):
    probabilities = []
    for i in range(1, VariantLength.max + 1):
        for j in range(1, VariantLength.max + 1):
            probabilities.append((
                (i, j),
                variant_length_probabilities_1[-i + VariantLength.index_offset] *
                variant_length_probabilities_2[j + VariantLength.index_offset] * extra_probability
            ))
    return probabilities


def maximum_variant_length_from(variant_length):
    if variant_length >= minimum_variant_length_that_need_infer:
        return maximum_variant_length_that_need_infer
    else:
        return variant_length


def quality_score_from(probability):
    """
    make a modification for quality score calculation. did not apply quality square for computation.
    """
    # p = probability
    # tmp = max(Phred_Trans * log(((1.0 - p) + 1e-300) / (p + 1e-300)) + 10, 0)
    # return float(round(tmp, 2))
    return probability

def possible_outcome_probabilites_with_indel_length_from(
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
        reference_base,
):
    homo_reference_probability = genotype_probabilities[Genotype.homo_reference]
    homo_variant_probability = genotype_probabilities[Genotype.homo_variant]
    hetero_variant_probability = genotype_probabilities[Genotype.hetero_variant]
    variant_length_0_probability = (
            variant_length_probabilities_1[0 + VariantLength.index_offset] *
            variant_length_probabilities_2[0 + VariantLength.index_offset]
    )

    reference_gt21 = gt21_enum_from_label(reference_base + reference_base)
    homo_Ref_probability = (
            variant_length_0_probability * homo_reference_probability * gt21_probabilities[reference_gt21]
    )

    homo_SNP_probabilities = [(
            variant_length_0_probability * homo_variant_probability * gt21_probabilities[gt21]
    ) for gt21 in HOMO_SNP_GT21]
    hetero_SNP_probabilities = [(
            variant_length_0_probability * hetero_variant_probability * gt21_probabilities[gt21]
    ) for gt21 in HETERO_SNP_GT21]

    # Insertion
    homo_Ins_lengths, homo_Ins_probabilities = zip(*homo_Ins_tuples_from(
        variant_length_probabilities_1, variant_length_probabilities_2,
        homo_variant_probability * gt21_probabilities[GT21_Type.InsIns]
    ))
    homo_Ins_lengths, homo_Ins_probabilities = list(homo_Ins_lengths), list(homo_Ins_probabilities)
    hetero_InsIns_length_tuples, hetero_InsIns_probabilities = zip(*hetero_InsIns_tuples_from(
        variant_length_probabilities_1, variant_length_probabilities_2,
        hetero_variant_probability * gt21_probabilities[GT21_Type.InsIns]
    ))
    hetero_InsIns_length_tuples, hetero_InsIns_probabilities = (
        list(hetero_InsIns_length_tuples), list(hetero_InsIns_probabilities)
    )
    hetero_ACGT_Ins_tuples = []
    gt21_base_tuples = [(GT21_Type.AIns, "A"), (GT21_Type.CIns, "C"), (GT21_Type.GIns, "G"), (GT21_Type.TIns, "T")]
    for length_tuples, p in hetero_Ins_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
        for gt21, hetero_base in gt21_base_tuples:
            hetero_ACGT_Ins_tuples.append((
                hetero_base,
                length_tuples,
                p * gt21_probabilities[gt21] * hetero_variant_probability
            ))
    hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities = zip(*hetero_ACGT_Ins_tuples)
    hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities = (
        list(hetero_ACGT_Ins_bases), list(hetero_ACGT_Ins_lengths), list(hetero_ACGT_Ins_probabilities)
    )

    # Deletion
    homo_Del_lengths, homo_Del_probabilities = zip(*homo_Del_tuples_from(
        variant_length_probabilities_1, variant_length_probabilities_2,
        homo_variant_probability * gt21_probabilities[GT21_Type.DelDel]
    ))
    homo_Del_lengths, homo_Del_probabilities = list(homo_Del_lengths), list(homo_Del_probabilities)
    hetero_DelDel_length_tuples, hetero_DelDel_probabilities = zip(*hetero_DelDel_tuples_from(
        variant_length_probabilities_1, variant_length_probabilities_2,
        hetero_variant_probability * gt21_probabilities[GT21_Type.DelDel]
    ))
    hetero_DelDel_length_tuples, hetero_DelDel_probabilities = (
        list(hetero_DelDel_length_tuples), list(hetero_DelDel_probabilities)
    )
    hetero_ACGT_Del_tuples = []
    gt21_base_tuples = [(GT21_Type.ADel, "A"), (GT21_Type.CDel, "C"), (GT21_Type.GDel, "G"), (GT21_Type.TDel, "T")]
    for length_tuples, p in hetero_Del_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
        for gt21, hetero_base in gt21_base_tuples:
            hetero_ACGT_Del_tuples.append((
                hetero_base,
                length_tuples,
                p * gt21_probabilities[gt21] * hetero_variant_probability
            ))
    hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities = zip(*hetero_ACGT_Del_tuples)
    hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities = (
        list(hetero_ACGT_Del_bases), list(hetero_ACGT_Del_lengths), list(hetero_ACGT_Del_probabilities)
    )

    # InsDel
    hetero_InsDel_length_tuples, hetero_InsDel_probabilities = zip(*hetero_InsDel_tuples_from(
        variant_length_probabilities_1, variant_length_probabilities_2,
        hetero_variant_probability * gt21_probabilities[GT21_Type.InsDel]
    ))
    hetero_InsDel_length_tuples, hetero_InsDel_probabilities = (
        list(hetero_InsDel_length_tuples), list(hetero_InsDel_probabilities)
    )

    return (
        homo_Ref_probability,
        homo_SNP_probabilities,
        hetero_SNP_probabilities,
        homo_Ins_lengths, homo_Ins_probabilities,
        hetero_InsIns_length_tuples, hetero_InsIns_probabilities,
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities,
        homo_Del_lengths, homo_Del_probabilities,
        hetero_DelDel_length_tuples, hetero_DelDel_probabilities,
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities,
        hetero_InsDel_length_tuples, hetero_InsDel_probabilities,
    )



def possible_outcome_probabilites_from(
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
        reference_base,
        alt_info_dict,
        add_indel_length=False,
):
    homo_reference_probability = genotype_probabilities[Genotype.homo_reference]
    homo_variant_probability = genotype_probabilities[Genotype.homo_variant]
    hetero_variant_probability = genotype_probabilities[Genotype.hetero_variant]

    reference_gt21 = gt21_enum_from_label(reference_base + reference_base)

    if not add_indel_length:
        homo_Ref_probability = (homo_reference_probability * gt21_probabilities[reference_gt21]
                                )
        homo_SNP_probabilities = [homo_variant_probability * gt21_probabilities[gt21]
                                  for gt21 in HOMO_SNP_GT21]
        hetero_SNP_probabilities = [hetero_variant_probability * gt21_probabilities[gt21]
                                    for gt21 in HETERO_SNP_GT21]
        if homo_reference_probability >= 0.5 and gt21_probabilities[
            reference_gt21] >= 0.5:
            return [homo_Ref_probability]
        # Insertion
        homo_Ins_probabilities = [homo_variant_probability * gt21_probabilities[GT21_Type.InsIns]]
        homo_Ins_lengths = []
        hetero_InsIns_probabilities = [hetero_variant_probability * gt21_probabilities[GT21_Type.InsIns]]
        hetero_InsIns_length_tuples = []
        hetero_ACGT_Ins_probabilities = []
        gt21_base_tuples = [(GT21_Type.AIns, "A"), (GT21_Type.CIns, "C"), (GT21_Type.GIns, "G"), (GT21_Type.TIns, "T")]
        for gt21, hetero_base in gt21_base_tuples:
            hetero_ACGT_Ins_probabilities.append(gt21_probabilities[gt21] * hetero_variant_probability)
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths = [], []

        # Deletion
        homo_Del_probabilities = [homo_variant_probability * gt21_probabilities[GT21_Type.DelDel]]
        homo_Del_lengths = []

        hetero_DelDel_probabilities = [hetero_variant_probability * gt21_probabilities[GT21_Type.DelDel]]
        hetero_DelDel_length_tuples = []

        hetero_ACGT_Del_probabilities = []
        gt21_base_tuples = [(GT21_Type.ADel, "A"), (GT21_Type.CDel, "C"), (GT21_Type.GDel, "G"), (GT21_Type.TDel, "T")]
        for gt21, hetero_base in gt21_base_tuples:
            hetero_ACGT_Del_probabilities.append(gt21_probabilities[gt21] * hetero_variant_probability
                                                 )
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths = [], []

        # InsDel
        hetero_InsDel_probabilities = [hetero_variant_probability * gt21_probabilities[GT21_Type.InsDel]]
        hetero_InsDel_length_tuples = []

    else:
        variant_length_0_probability_1 = variant_length_probabilities_1[0 + VariantLength.index_offset]
        variant_length_0_probability_2 = variant_length_probabilities_2[0 + VariantLength.index_offset]
        variant_length_0_probability = (variant_length_0_probability_1 * variant_length_0_probability_2)

        reference_gt21 = gt21_enum_from_label(reference_base + reference_base)
        homo_Ref_probability = (
                variant_length_0_probability * homo_reference_probability * gt21_probabilities[reference_gt21]
        )
        if variant_length_0_probability_1 >= 0.5 and variant_length_0_probability_2 >= 0.5 and homo_reference_probability >= 0.5 and \
                gt21_probabilities[
                    reference_gt21] >= 0.5:
            return [homo_Ref_probability]

        homo_SNP_probabilities = [(
                variant_length_0_probability * homo_variant_probability * gt21_probabilities[gt21]
        ) for gt21 in HOMO_SNP_GT21]
        hetero_SNP_probabilities = [(
                variant_length_0_probability * hetero_variant_probability * gt21_probabilities[gt21]
        ) for gt21 in HETERO_SNP_GT21]

        # Insertion
        homo_Ins_lengths, homo_Ins_probabilities = zip(*homo_Ins_tuples_from(
            variant_length_probabilities_1, variant_length_probabilities_2,
            homo_variant_probability * gt21_probabilities[GT21_Type.InsIns]
        ))
        homo_Ins_lengths, homo_Ins_probabilities = list(homo_Ins_lengths), list(homo_Ins_probabilities)
        hetero_InsIns_length_tuples, hetero_InsIns_probabilities = zip(*hetero_InsIns_tuples_from(
            variant_length_probabilities_1, variant_length_probabilities_2,
            hetero_variant_probability * gt21_probabilities[GT21_Type.InsIns]
        ))
        hetero_InsIns_length_tuples, hetero_InsIns_probabilities = (
            list(hetero_InsIns_length_tuples), list(hetero_InsIns_probabilities)
        )
        hetero_ACGT_Ins_tuples = []
        gt21_base_tuples = [(GT21_Type.AIns, "A"), (GT21_Type.CIns, "C"), (GT21_Type.GIns, "G"), (GT21_Type.TIns, "T")]
        for length_tuples, p in hetero_Ins_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
            for gt21, hetero_base in gt21_base_tuples:
                hetero_ACGT_Ins_tuples.append((
                    hetero_base,
                    length_tuples,
                    p * gt21_probabilities[gt21] * hetero_variant_probability
                ))
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities = zip(*hetero_ACGT_Ins_tuples)
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities = (
            list(hetero_ACGT_Ins_bases), list(hetero_ACGT_Ins_lengths), list(hetero_ACGT_Ins_probabilities)
        )

        # Deletion
        homo_Del_lengths, homo_Del_probabilities = zip(*homo_Del_tuples_from(
            variant_length_probabilities_1, variant_length_probabilities_2,
            homo_variant_probability * gt21_probabilities[GT21_Type.DelDel]
        ))
        homo_Del_lengths, homo_Del_probabilities = list(homo_Del_lengths), list(homo_Del_probabilities)
        hetero_DelDel_length_tuples, hetero_DelDel_probabilities = zip(*hetero_DelDel_tuples_from(
            variant_length_probabilities_1, variant_length_probabilities_2,
            hetero_variant_probability * gt21_probabilities[GT21_Type.DelDel]
        ))
        hetero_DelDel_length_tuples, hetero_DelDel_probabilities = (
            list(hetero_DelDel_length_tuples), list(hetero_DelDel_probabilities)
        )
        hetero_ACGT_Del_tuples = []
        gt21_base_tuples = [(GT21_Type.ADel, "A"), (GT21_Type.CDel, "C"), (GT21_Type.GDel, "G"), (GT21_Type.TDel, "T")]
        for length_tuples, p in hetero_Del_tuples_from(variant_length_probabilities_1, variant_length_probabilities_2):
            for gt21, hetero_base in gt21_base_tuples:
                hetero_ACGT_Del_tuples.append((
                    hetero_base,
                    length_tuples,
                    p * gt21_probabilities[gt21] * hetero_variant_probability
                ))
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities = zip(*hetero_ACGT_Del_tuples)
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities = (
            list(hetero_ACGT_Del_bases), list(hetero_ACGT_Del_lengths), list(hetero_ACGT_Del_probabilities)
        )

        # InsDel
        hetero_InsDel_length_tuples, hetero_InsDel_probabilities = zip(*hetero_InsDel_tuples_from(
            variant_length_probabilities_1, variant_length_probabilities_2,
            hetero_variant_probability * gt21_probabilities[GT21_Type.InsDel]
        ))
        hetero_InsDel_length_tuples, hetero_InsDel_probabilities = (
            list(hetero_InsDel_length_tuples), list(hetero_InsDel_probabilities)
        )

    return (
        homo_Ref_probability,
        homo_SNP_probabilities,
        hetero_SNP_probabilities,
        homo_Ins_lengths, homo_Ins_probabilities,
        hetero_InsIns_length_tuples, hetero_InsIns_probabilities,
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities,
        homo_Del_lengths, homo_Del_probabilities,
        hetero_DelDel_length_tuples, hetero_DelDel_probabilities,
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities,
        hetero_InsDel_length_tuples, hetero_InsDel_probabilities,
    )


def argmax(l):
    return max(zip(l, range(len(l))))[1]
def find_alt_base(alt_info_dict, alternate_base=None):
    # double check whether alternate base exists, depth gap may happen when candidate depth is extreme high to infer.
    max_depth_gap = 9
    sorted_alt_bases = sorted([(alt_base[1], count) for alt_base, count in alt_info_dict.items() if alt_base[0] == 'X'],
                              key=lambda x: x[1], reverse=True)
    alt_count = [item[1] for item in sorted_alt_bases if item[0] == alternate_base]
    if not len(sorted_alt_bases):
        return [], None
    if not len(alt_count) or sorted_alt_bases[0][1] - alt_count[0] >= max_depth_gap:
        alternate_base = sorted_alt_bases[0][0]  #
    sorted_alt_bases = [item[0] for item in sorted_alt_bases]
    return sorted_alt_bases, alternate_base


def output_from(
        reference_sequence,
        contig,
        position,
        tensor_position_center,
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
        output_config,
        output_utilities,
        alt_info_dict
):
    add_indel_length = output_config.add_indel_length
    reference_base_ACGT = BASE2ACGT[reference_sequence[tensor_position_center]]

    all_pro = possible_outcome_probabilites_from(
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
        reference_base=reference_base_ACGT,
        alt_info_dict=alt_info_dict,
        add_indel_length=add_indel_length,
    )

    # return if a candidate is homo reference
    if len(all_pro) == 1:
        return (
            (True, False, False, False, False, False, False, False, False, False),
            (reference_base_ACGT, reference_base_ACGT), (all_pro[0])
        )
    (
        homo_Ref_probability,
        homo_SNP_probabilities,
        hetero_SNP_probabilities,
        homo_Ins_lengths, homo_Ins_probabilities,
        hetero_InsIns_length_tuples, hetero_InsIns_probabilities,
        hetero_ACGT_Ins_bases, hetero_ACGT_Ins_lengths, hetero_ACGT_Ins_probabilities,
        homo_Del_lengths, homo_Del_probabilities,
        hetero_DelDel_length_tuples, hetero_DelDel_probabilities,
        hetero_ACGT_Del_bases, hetero_ACGT_Del_lengths, hetero_ACGT_Del_probabilities,
        hetero_InsDel_length_tuples, hetero_InsDel_probabilities,
    ) = all_pro
    maximum_probability = 0.0
    maximum_loops = 4
    reference_base, alternate_base = None, None
    loop_index = 0
    while (reference_base is None or alternate_base is None) and loop_index < maximum_loops:
        loop_index += 1
        maximum_probability = max(
            homo_Ref_probability,
            max(homo_SNP_probabilities),
            max(hetero_SNP_probabilities),
            max(homo_Ins_probabilities) if len(homo_Ins_probabilities) else 0,
            max(homo_Del_probabilities) if len(homo_Del_probabilities) else 0,
            max(hetero_ACGT_Ins_probabilities) if len(hetero_ACGT_Ins_probabilities) else 0,
            max(hetero_InsIns_probabilities) if len(hetero_InsIns_probabilities) else 0,
            max(hetero_ACGT_Del_probabilities) if len(hetero_ACGT_Del_probabilities) else 0,
            max(hetero_DelDel_probabilities) if len(hetero_DelDel_probabilities) else 0,
            max(hetero_InsDel_probabilities) if len(hetero_InsDel_probabilities) else 0,
        )

        is_reference = maximum_probability == homo_Ref_probability
        if is_reference:
            return (
                (True, False, False, False, False, False, False, False, False, False),
                (reference_base_ACGT, reference_base_ACGT), (maximum_probability)
            )

        is_homo_SNP = maximum_probability in homo_SNP_probabilities
        is_hetero_SNP = maximum_probability in hetero_SNP_probabilities
        is_homo_insertion = maximum_probability in homo_Ins_probabilities
        is_hetero_ACGT_Ins = maximum_probability in hetero_ACGT_Ins_probabilities
        is_hetero_InsIns = maximum_probability in hetero_InsIns_probabilities
        is_homo_deletion = maximum_probability in homo_Del_probabilities
        is_hetero_ACGT_Del = maximum_probability in hetero_ACGT_Del_probabilities
        is_hetero_DelDel = maximum_probability in hetero_DelDel_probabilities
        is_insertion_and_deletion = maximum_probability in hetero_InsDel_probabilities

        if is_homo_SNP:
            reference_base = reference_sequence[tensor_position_center]
            base1, base2 = homo_SNP_bases_from(homo_SNP_probabilities)
            alternate_base = base1 if base1 != reference_base else base2
            sorted_alt_bases, alternate_base = find_alt_base(alt_info_dict, alternate_base)

        elif is_hetero_SNP:
            base1, base2 = hetero_SNP_bases_from(hetero_SNP_probabilities)
            reference_base = reference_sequence[tensor_position_center]
            is_multi = base1 != reference_base and base2 != reference_base
            if is_multi:
                sorted_alt_bases, _ = find_alt_base(alt_info_dict)
                if len(sorted_alt_bases) == 0:
                    break
                if len(sorted_alt_bases) < 2:
                    alternate_base = sorted_alt_bases[0]
                    hetero_SNP_probabilities[np.argmax(hetero_SNP_probabilities)] = 0.0
                    break
                alternate_base = ','.join(sorted_alt_bases[:2])
            else:
                alternate_base = base1 if base1 != reference_base else base2
                sorted_alt_bases, alternate_base = find_alt_base(alt_info_dict, alternate_base)


        elif is_homo_insertion:
            variant_length = None
            if add_indel_length:
                idx = homo_Ins_probabilities.index(maximum_probability)
                variant_length = homo_Ins_lengths[idx]
            insertion_bases = insertion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_insertion_length=variant_length if variant_length and variant_length < VariantLength.max else None,
            )

            insertion_length = len(insertion_bases)
            if insertion_length == 0:
                break
            reference_base = reference_sequence[tensor_position_center]
            alternate_base = insertion_bases

        elif is_hetero_ACGT_Ins:
            idx = hetero_ACGT_Ins_probabilities.index(maximum_probability)
            variant_length = None
            if add_indel_length:
                hetero_Ins_base = hetero_ACGT_Ins_bases[idx]
                variant_length = hetero_ACGT_Ins_lengths[idx]
            else:
                hetero_Ins_base = ACGT[idx]
            insertion_bases = insertion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_insertion_length=variant_length if variant_length and variant_length < VariantLength.max else None,

            )
            insertion_length = len(insertion_bases)
            if insertion_length == 0:
                break
            reference_base = reference_sequence[tensor_position_center]
            alternate_base = insertion_bases

            is_SNP_Ins_multi = hetero_Ins_base != reference_base
            if is_SNP_Ins_multi:
                sorted_alt_bases, _ = find_alt_base(alt_info_dict)
                if len(sorted_alt_bases) == 0:
                    break
                else:
                    alternate_base = "{},{}".format(sorted_alt_bases[0], alternate_base)

        elif is_hetero_InsIns:
            insertion_bases_list = []
            if add_indel_length:
                idx = hetero_InsIns_probabilities.index(maximum_probability)
                variant_length_1, variant_length_2 = hetero_InsIns_length_tuples[idx]
                del hetero_InsIns_probabilities[idx]
                del hetero_InsIns_length_tuples[idx]

                insertion_bases1 = insertion_bases_using_alt_info_from(
                    alt_info_dict=alt_info_dict,
                    propose_insertion_length=variant_length_1 if variant_length_1 and variant_length_1 < VariantLength.max else None,
                )
                if len(insertion_bases1):
                    insertion_bases2 = insertion_bases_using_alt_info_from(
                        alt_info_dict=alt_info_dict,
                        propose_insertion_length=variant_length_2 if variant_length_2 and variant_length_2 < VariantLength.max else None,
                        insertion_bases_to_ignore=insertion_bases1
                    )
                    if len(insertion_bases2):
                        insertion_bases_list = [insertion_bases1, insertion_bases2]
                if len(insertion_bases_list) < 2:
                    insertion_bases_list = insertion_bases_using_alt_info_from(
                        alt_info_dict=alt_info_dict,
                        return_multi=True
                    )
            else:
                insertion_bases_list = insertion_bases_using_alt_info_from(
                    alt_info_dict=alt_info_dict,
                    return_multi=True
                )
            if len(insertion_bases_list) < 2:
                break
            insertion_bases, another_insertion_bases = insertion_bases_list

            reference_base = reference_sequence[tensor_position_center]
            alternate_base = insertion_bases

            alternate_base_1 = another_insertion_bases
            alternate_base_2 = alternate_base
            if alternate_base_1 != alternate_base_2:
                alternate_base = "{},{}".format(alternate_base_1, alternate_base_2)
            else:
                reference_base, alternate_base = None, None

        elif is_homo_deletion:
            variant_length = None
            if add_indel_length:
                idx = homo_Del_probabilities.index(maximum_probability)
                variant_length = homo_Del_lengths[idx]

            deletion_bases = deletion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_deletion_length=variant_length if variant_length and variant_length < VariantLength.max else None,
            )
            deletion_length = len(deletion_bases)
            if deletion_length == 0:
                break
            reference_base = reference_sequence[tensor_position_center] + deletion_bases
            alternate_base = reference_base[0]

        elif is_hetero_ACGT_Del:
            variant_length = None
            idx = hetero_ACGT_Del_probabilities.index(maximum_probability)
            if add_indel_length:
                variant_length = hetero_ACGT_Del_lengths[idx]
                hetero_Del_base = hetero_ACGT_Del_bases[idx]
            else:
                hetero_Del_base = ACGT[idx]
            deletion_bases = deletion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_deletion_length=variant_length if variant_length and variant_length < VariantLength.max else None,
            )
            deletion_length = len(deletion_bases)
            if deletion_length == 0:
                break
            reference_base = reference_sequence[tensor_position_center] + deletion_bases
            alternate_base = reference_base[0]

            is_SNP_Del_multi = hetero_Del_base != reference_base[0]
            if is_SNP_Del_multi:
                alternate_base_1 = alternate_base
                alternate_base_2 = hetero_Del_base + reference_base[1:]
                alternate_base = "{},{}".format(alternate_base_1, alternate_base_2)

        elif is_hetero_DelDel:
            deletion_bases_list = []
            if add_indel_length:
                idx = hetero_DelDel_probabilities.index(maximum_probability)
                variant_length_1, variant_length_2 = sorted(hetero_DelDel_length_tuples[idx],
                                                            reverse=True)  # longer deletion should be in first position
                deletion_base1 = deletion_bases_using_alt_info_from(
                    alt_info_dict=alt_info_dict,
                    propose_deletion_length=variant_length_1 if variant_length_1 and variant_length_1 < VariantLength.max else None,
                )
                if len(deletion_base1) > 0:
                    deletion_base2 = deletion_bases_using_alt_info_from(
                        alt_info_dict=alt_info_dict,
                        propose_deletion_length=variant_length_2 if variant_length_2 and variant_length_2 < VariantLength.max else None,
                        deletion_bases_to_ignore=deletion_base1
                    )
                    if len(deletion_base2) > 0:
                        deletion_bases_list = [deletion_base1, deletion_base2] if len(deletion_base1) > len(
                            deletion_base2) else [deletion_base2, deletion_base1]
                if len(deletion_bases_list) < 2:
                    deletion_bases_list = deletion_bases_using_alt_info_from(
                        return_multi=True,
                        alt_info_dict=alt_info_dict
                    )
            else:
                deletion_bases_list = deletion_bases_using_alt_info_from(
                    return_multi=True,
                    alt_info_dict=alt_info_dict
                )

            if len(deletion_bases_list) < 2:
                break

            deletion_bases, deletion_bases1 = deletion_bases_list

            reference_base = reference_sequence[tensor_position_center] + deletion_bases
            alternate_base = reference_base[0]

            alternate_base_1 = alternate_base
            alternate_base_2 = reference_base[0] + reference_base[len(deletion_bases1) + 1:]
            if (
                    alternate_base_1 != alternate_base_2 and
                    reference_base != alternate_base_1 and reference_base != alternate_base_2
            ):
                alternate_base = "{},{}".format(alternate_base_1, alternate_base_2)
            else:
                reference_base, alternate_base = None, None

        elif is_insertion_and_deletion:
            variant_length_1, variant_length_2 = None, None
            if add_indel_length:
                idx = hetero_InsDel_probabilities.index(maximum_probability)
                variant_length_1, variant_length_2 = hetero_InsDel_length_tuples[idx]

            insertion_bases = insertion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_insertion_length=variant_length_2 if variant_length_2 and variant_length_2 < VariantLength.max else None
            )
            insertion_length = len(insertion_bases)

            deletion_bases = deletion_bases_using_alt_info_from(
                alt_info_dict=alt_info_dict,
                propose_deletion_length=variant_length_1 if variant_length_1 and variant_length_1 < VariantLength.max else None,
            )
            deletion_length = len(deletion_bases)

            if insertion_length == 0 or deletion_length == 0:
                break
            reference_base = reference_sequence[tensor_position_center] + deletion_bases
            alternate_base = "{},{}".format(
                reference_base[0],
                insertion_bases + reference_base[1:]
            )

    return (
               is_reference, is_homo_SNP, is_hetero_SNP,
               is_homo_insertion, is_hetero_ACGT_Ins, is_hetero_InsIns,
               is_homo_deletion, is_hetero_ACGT_Del, is_hetero_DelDel,
               is_insertion_and_deletion
           ), (reference_base, alternate_base), (maximum_probability)


def batch_output_for_ensemble(X, batch_chr_pos_seq, alt_info_list, batch_Y, output_config, output_utilities):
    batch_size = len(batch_chr_pos_seq)
    batch_gt21_probabilities, batch_genotype_probabilities, = batch_Y

    if len(batch_gt21_probabilities) != batch_size:
        sys.exit(
            "Inconsistent shape between input tensor and output predictions %d/%d" %
            (batch_size, len(batch_gt21_probabilities))
        )

    tensor_position_center = param.flankingBaseNum

    for (
            x,
            chr_pos_seq,
            gt21_probabilities,
            genotype_probabilities,
            alt_info
    ) in zip(
        X,
        batch_chr_pos_seq,
        batch_gt21_probabilities,
        batch_genotype_probabilities,
        alt_info_list
    ):
        if output_config.tensor_fn != 'PIPE':
            chromosome, position, reference_sequence = chr_pos_seq.decode().rstrip().split(":")
        else:
            chromosome, position, reference_sequence = chr_pos_seq

        position = int(position)

        if reference_sequence[tensor_position_center] not in BASIC_BASES:
            continue

        output_utilities.output(
            "\t".join(
                [
                    chromosome,
                    str(position),
                    reference_sequence,
                    alt_info.decode(),
                    ' '.join(["{:0.6f}".format(p) for p in list(gt21_probabilities)]),
                    ' '.join(["{:0.6f}".format(p) for p in list(genotype_probabilities)])]
            )
        )

def output_vcf_from_probability(
        chromosome,
        position,
        reference_base,
        normal_alt_info,
        tumor_alt_info,
        gt21_probabilities,
        genotype_probabilities=None,
        variant_length_probabilities_1=None,
        variant_length_probabilities_2=None,
        output_config=None,
        vcf_writer=None,
):
    # if type(chr_pos_seq) == np.memmap:
    #     chr_pos_seq = chr_pos_seq[0].decode()
    # elif type(chr_pos_seq) == np.bytes_ or type(chr_pos_seq) == bytes:
    #     chr_pos_seq = chr_pos_seq.decode()
    #
    # chromosome, position, reference_sequence = chr_pos_seq.rstrip().split(':')
    # position = int(position)
    #
    # tensor_position_center = param.flankingBaseNum
    # information_string = "P" if output_config.pileup else 'F'
    #
    # if type(alt_info) == np.memmap:
    #     alt_info = alt_info[0].decode()
    # elif type(alt_info) == np.bytes_ or type(alt_info) == bytes:
    #     alt_info = alt_info.decode()
    def decode_alt_info(alt_info):
        alt_info = alt_info.rstrip().split('-')
        read_depth = int(alt_info[0])  # alt_info
        indel_str = alt_info[1] if len(alt_info) > 1 else ''
        seqs = indel_str.split(' ')
        alt_info_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
        return alt_info_dict, read_depth

    normal_alt_info_dict, normal_read_depth = decode_alt_info(normal_alt_info)
    tumor_alt_info_dict, tumor_read_depth = decode_alt_info(tumor_alt_info)

    # output_info = output_from(
    #     reference_sequence,
    #     chromosome,
    #     position,
    #     tensor_position_center,
    #     gt21_probabilities,
    #     genotype_probabilities,
    #     variant_length_probabilities_1,
    #     variant_length_probabilities_2,
    #     output_config,
    #     output_utilities,
    #     alt_info_dict
    # )

    somatic_arg_index = param.somatic_arg_index
    alternate_base = reference_base
    arg_index = argmax(gt21_probabilities)
    is_reference = arg_index == 0
    is_germline = arg_index == 1
    is_tumor = arg_index == somatic_arg_index
    maximum_probability = gt21_probabilities[arg_index]
    def rank_somatic_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth):
        support_alt_dict = {}
        for tumor_alt, tumor_count in tumor_alt_info_dict.items():
            tumor_af = tumor_count / float(tumor_read_depth)
            normal_count = normal_alt_info_dict[tumor_alt] if tumor_alt in normal_alt_info_dict else 0
            normal_af = normal_count / float(normal_read_depth) if normal_read_depth > 0 else 0
            if tumor_af - normal_af > 0:
                support_alt_dict[tumor_alt] = tumor_af - normal_af
        if len(support_alt_dict) == 0:
            return "", 0, 0
        alt_type_list = sorted(support_alt_dict.items(), key=lambda x: x[1], reverse=True)
        best_match_alt = alt_type_list[0][0]
        tumor_supported_reads_count = tumor_alt_info_dict[best_match_alt]
        normal_supported_reads_count = normal_alt_info_dict[best_match_alt] if best_match_alt in normal_alt_info_dict else 0
        return best_match_alt, tumor_supported_reads_count, normal_supported_reads_count

    def rank_germline_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth):
        support_alt_dict = {}
        for tumor_alt, tumor_count in tumor_alt_info_dict.items():
            tumor_af = tumor_count / float(tumor_read_depth)
            normal_count = normal_alt_info_dict[tumor_alt] if tumor_alt in normal_alt_info_dict else 0
            normal_af = normal_count / float(normal_read_depth) if normal_read_depth > 0 else 0
            support_alt_dict[tumor_alt] = tumor_af + normal_af
        if len(support_alt_dict) == 0:
            return "", 0, 0
        alt_type_list = sorted(support_alt_dict.items(), key=lambda x: x[1], reverse=True)
        best_match_alt = alt_type_list[0][0]
        tumor_supported_reads_count = tumor_alt_info_dict[best_match_alt]
        normal_supported_reads_count = normal_alt_info_dict[
            best_match_alt] if best_match_alt in normal_alt_info_dict else 0
        return best_match_alt, tumor_supported_reads_count, normal_supported_reads_count


    if is_tumor:
        if tumor_read_depth <= 0:
            print("low tumor coverage")
            return
        best_match_alt, tumor_supported_reads_count, normal_supported_reads_count = rank_somatic_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth)

        if best_match_alt == "":
            return
        if best_match_alt[0] == 'X':
            alternate_base = best_match_alt[1]
            is_SNP = True
        elif best_match_alt[0] == 'I':
            alternate_base = best_match_alt[1:]
            is_INS = True
        elif best_match_alt[0] == 'D':
            alternate_base = reference_base
            reference_base = best_match_alt[1:]

    if is_germline and output_config.is_show_germline:
        best_match_alt, tumor_supported_reads_count, normal_supported_reads_count = rank_germline_alt(tumor_alt_info_dict, normal_alt_info_dict, tumor_read_depth, normal_read_depth)
        if best_match_alt == "":
            return
        if best_match_alt[0] == 'X':
            alternate_base = best_match_alt[1]
            is_SNP = True
        elif best_match_alt[0] == 'I':
            alternate_base = best_match_alt[1:]
            is_INS = True
        elif best_match_alt[0] == 'D':
            alternate_base = reference_base
            reference_base = best_match_alt[1:]

    if (not output_config.is_show_reference and is_reference) or (not is_reference and reference_base == alternate_base):
        return

    if (not output_config.is_show_germline and is_germline):
        return

    if reference_base is None or alternate_base is None:
        return

    # genotype string
    if is_reference:
        genotype_string = '0/0' #genotype_string_from(Genotype.homo_reference)

    else:
        genotype_string = "1/1"
        # allele frequency / supported reads

    def decode_alt_info(alt_info_dict, read_depth):
        alt_type_list = [{}, {}, {}]  # SNP I D
        snp_num, ins_num, del_num = 0, 0, 0
        for alt_type, count in alt_info_dict.items():
            count = int(count)
            if alt_type[0] == 'X':
                alt_type_list[0][alt_type[1]] = count
                snp_num += count
            elif alt_type[0] == 'I':
                alt_type_list[1][alt_type[1:]] = count
                ins_num += count
            elif alt_type[0] == 'D':
                alt_type_list[2][alt_type[1:]] = count
                del_num += count
        ref_num = max(0, read_depth - snp_num - ins_num - del_num)
        return alt_type_list, ref_num, snp_num, ins_num, del_num


    normal_alt_type_list, normal_ref_num, normal_snp_num, normal_ins_num, normal_del_num = decode_alt_info(alt_info_dict=normal_alt_info_dict, read_depth=normal_read_depth)
    tumor_alt_type_list, tumor_ref_num, tumor_snp_num, tumor_ins_num, tumor_del_num = decode_alt_info(alt_info_dict=tumor_alt_info_dict, read_depth=tumor_read_depth)

    if is_reference:
        normal_supported_reads_count = normal_ref_num
        tumor_supported_reads_count = tumor_ref_num
        alternate_base = "."

    normal_allele_frequency = min((normal_supported_reads_count / normal_read_depth) if normal_read_depth != 0 else 0.0, 1.0)
    tumor_allele_frequency = min((tumor_supported_reads_count / tumor_read_depth) if tumor_read_depth != 0 else 0.0, 1.0)
    

    # quality score
    quality_score = quality_score_from(maximum_probability)

    # filtration value
    filtration_value = filtration_value_from(
        quality_score_for_pass=output_config.quality_score_for_pass,
        quality_score=quality_score,
        is_reference=is_reference,
        is_germline=is_germline
    )

    information_string = "."

    read_depth = normal_read_depth + tumor_read_depth

    vcf_writer.write_row(CHROM=chromosome,
                         POS=position,
                         REF=reference_base,
                         ALT=alternate_base,
                         QUAL=quality_score,
                         FILTER=filtration_value,
                         INFO=information_string,
                         GT=genotype_string,
                         DP=read_depth,
                         NDP=normal_read_depth,
                         TDP=tumor_read_depth,
                         AF=tumor_allele_frequency,
                         NAF=normal_allele_frequency)
        #
        #
        # "%s\t%s\t.\t%s\t%s\t%.2f\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
        # chromosome,
        # position,
        # reference_base,
        # alternate_base,
        # quality_score,
        # filtration_value,
        # information_string,
        # genotype_string,
        # quality_score,
        # tumor_read_depth,
        # allele_frequency
    # ))

def call_variants_from_probability(args):
    # args.predict_fn = "/mnt/bal36/zxzheng/TMP/tmp_alt"
    # args.call_fn ="/mnt/bal36/zxzheng/TMP/tmp.vcf"
    output_config = OutputConfig(
        is_show_reference=args.show_ref,
        is_show_germline=args.show_germline,
        is_haploid_precise_mode_enabled=args.haploid_precise,
        is_haploid_sensitive_mode_enabled=args.haploid_sensitive,
        is_output_for_ensemble=args.output_for_ensemble,
        quality_score_for_pass=args.qual,
        tensor_fn=args.tensor_fn,
        input_probabilities=args.input_probabilities,
        add_indel_length=args.add_indel_length,
        gvcf=args.gvcf,
        pileup=args.pileup
    )

    call_fn = args.call_fn
    if call_fn != "PIPE":
        call_dir = os.path.dirname(call_fn)
        if not os.path.exists(call_dir):
            output = run("mkdir -p {}".format(call_dir), shell=True)
    vcf_writer = VcfWriter(vcf_fn=args.call_fn,
                           ref_fn=args.ref_fn,
                           ctg_name=args.ctg_name,
                           show_ref_calls=args.show_ref,
                           sample_name=args.sample_name,
                           )

    logging.info("Calling variants ...")
    variant_call_start_time = time()
    prediction_path = args.predict_fn
    # batch_output_method = batch_output


    # prediction_path = args.tensor_fn + '.prediction'
    if prediction_path != "PIPE":
        if not os.path.exists(prediction_path):
            print("[ERROR] Prediction path not found!")
            return
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, prediction_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    # output_utilities.gen_output_file()
    # output_utilities.output_header()
    for row_id, row in enumerate(fo):
        row = row.rstrip().split('\t')
        chromosome, position, reference_base, normal_alt_info, tumor_alt_info, prediction = row[:6]
        gt21_probabilities = [float(item) for item in prediction.split()]
        output_vcf_from_probability(
            chromosome,
            position,
            reference_base,
            normal_alt_info,
            tumor_alt_info,
            gt21_probabilities,
            output_config=output_config,
            vcf_writer=vcf_writer)

    logging.info("Total time elapsed: %.2f s" % (time() - variant_call_start_time))

    vcf_writer.close()
    # remove file if on variant in output
    if os.path.exists(args.call_fn):
        vcf_file = open(args.call_fn, 'r').readlines()
        if not len(vcf_file):
            os.remove(args.call_fn)
        for row in vcf_file:
            if row[0] != '#':
                return
        logging.info("[INFO] No vcf output for file {}, remove empty file".format(args.call_fn))
        os.remove(args.call_fn)


def DataGenerator(dataset, num_epoch, batch_size, chunk_start_pos, chunk_end_pos):
    for idx in range(num_epoch):
        start_pos = chunk_start_pos + idx * batch_size
        end_pos = min(chunk_start_pos + (idx + 1) * batch_size, chunk_end_pos)
        position_matrix = dataset.position_matrix[start_pos:end_pos]
        position = dataset.position[start_pos:end_pos]  # .flatten()
        alt_info_list = dataset.alt_info[start_pos:end_pos]  # .flatten()
        yield position_matrix, position, alt_info_list


def main():
    parser = ArgumentParser(description="Call variants using a trained model and tensors of candidate variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a trained model for variant calling, required")

    parser.add_argument('--call_fn', type=str, default="clair3",
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required if --gvcf is enabled")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--ctg_start', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctg_end', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=0,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    # options for advanced users
    parser.add_argument('--temp_file_dir', type=str, default='./',
                        help="EXPERIMENTAL: The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--haploid_precise', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant")

    parser.add_argument('--haploid_sensitive', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant")

    # options for debug purpose
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building your own pipeline")

    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    parser.add_argument('--input_probabilities', action='store_true',
                        help="DEBUG: Use network probability outputs as input and generate variants from them")

    parser.add_argument('--output_probabilities', action='store_true',
                        help="DEBUG: Output the network probabilities of gt21, genotype, indel_length_1 and indel_length_2")

    # options for internal process control
    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Include indel length in training and calling, false for pileup and true for raw alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Generating outputs for ensemble model calling
    parser.add_argument('--output_for_ensemble', action='store_true',
                        help=SUPPRESS)

    ## Use bin file from pytables to speed up calling.
    parser.add_argument('--is_from_tables', action='store_true',
                        help=SUPPRESS)

    ## Output reference calls
    parser.add_argument('--show_ref', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--show_germline', action='store_true',
                        help=SUPPRESS)
    args = parser.parse_args()

    # if len(sys.argv[1:]) == 0:
    #     parser.print_help()
    #     sys.exit(1)

    call_variants_from_probability(args)


if __name__ == "__main__":
    main()
