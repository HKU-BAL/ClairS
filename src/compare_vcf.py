import os
import subprocess
import shlex
from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict


from shared.utils import log_error, log_warning, file_path_from, subprocess_popen
from shared.vcf import VcfReader, VcfWriter
from shared.interval_tree import bed_tree_from, is_region_in
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def cal_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return round(precision, 6), round(recall, 6), round(f1_score, 6)


def compare_vcf(args):
    """
    How som.py works
    """
    output_fn = args.output_fn
    output_dir = args.output_dir
    truth_vcf_fn = args.truth_vcf_fn
    input_vcf_fn = args.input_vcf_fn
    bed_fn = args.bed_fn
    sample_name = args.sampleName
    ref_fn = args.ref_fn
    ctg_name = args.ctg_name
    skip_genotyping = args.skip_genotyping
    input_filter_tag = args.input_filter_tag
    truth_filter_tag = args.truth_filter_tag
    fp_bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=ctg_name)
    # fp_bed_tree = {}

    truth_vcf_reader = VcfReader(vcf_fn=truth_vcf_fn, ctg_name=ctg_name, show_ref=False,keep_row_str=True, filter_tag=truth_filter_tag)
    truth_vcf_reader.read_vcf()
    truth_variant_dict = truth_vcf_reader.variant_dict

    input_vcf_reader = VcfReader(vcf_fn=input_vcf_fn, ctg_name=ctg_name, show_ref=False,keep_row_str=True, filter_tag=input_filter_tag)
    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict

    if output_fn:
        output_file = open(output_fn, 'w')
    else:
        output_file = None

    tp_snp, tp_ins, tp_del, fp_snp, fp_ins, fp_del, fn_snp, fn_ins, fn_del, fp_snp_truth, fp_ins_truth, fp_del_truth = 0,0,0,0,0,0,0,0,0,0,0,0
    truth_set = set()
    truth_snp, truth_ins, truth_del = 0,0,0
    query_snp, query_ins, query_del = 0,0,0
    pos_out_of_bed = 0

    fp_set = set()
    fn_set = set()
    fp_fn_set = set()
    tp_set = set()
    for key, vcf_infos in input_variant_dict.items():
        pos = key if ctg_name is not None else key[1]
        contig = ctg_name if ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                    contig_name=contig,
                                                    region_start=pos-1,
                                                    region_end=pos)
        if not pass_bed_region:
            pos_out_of_bed += 1
            # print(pos)
            continue

        ref_base = vcf_infos.reference_bases
        alt_base = vcf_infos.alternate_bases[0]
        genotype = vcf_infos.genotype
        qual = vcf_infos.qual
        is_snp = len(ref_base) == 1 and len(alt_base) == 1
        is_ins = len(ref_base) < len(alt_base)
        is_del = len(ref_base) > len(alt_base)

        if key not in truth_variant_dict and genotype != (0, 0):
            fp_snp = fp_snp + 1 if is_snp else fp_snp
            fp_ins = fp_ins + 1 if is_ins else fp_ins
            fp_del = fp_del + 1 if is_del else fp_del
            if fp_snp:
                fp_set.add(pos)

        if key in truth_variant_dict:
            vcf_infos = truth_variant_dict[key]
            truth_ref_base = vcf_infos.reference_bases
            truth_alt_base = vcf_infos.alternate_bases[0]
            truth_genotype = vcf_infos.genotype
            is_snp_truth = len(truth_ref_base) == 1 and len(truth_alt_base) == 1
            is_ins_truth = len(truth_ref_base) < len(truth_alt_base)
            is_del_truth = len(truth_ref_base) > len(truth_alt_base)

            if genotype == (0, 0) and truth_genotype == (0, 0):
                continue

            genotype_match = skip_genotyping or (truth_genotype == genotype)
            if truth_ref_base == ref_base and truth_alt_base == alt_base and genotype_match:
                tp_snp = tp_snp + 1 if is_snp else tp_snp
                tp_ins = tp_ins + 1 if is_ins else tp_ins
                tp_del = tp_del + 1 if is_del else tp_del
                if tp_snp or is_snp_truth:
                    tp_set.add(pos)
            else:
                fp_snp = fp_snp + 1 if is_snp else fp_snp
                fp_ins = fp_ins + 1 if is_ins else fp_ins
                fp_del = fp_del + 1 if is_del else fp_del

                fn_snp = fn_snp + 1 if is_snp_truth else fn_snp
                fn_ins = fn_ins + 1 if is_ins_truth else fn_ins
                fn_del = fn_del + 1 if is_del_truth else fn_del

                if fn_snp or fp_snp:
                    fp_fn_set.add(pos)

            truth_set.add(pos)

    for key, vcf_infos in truth_variant_dict.items():
        pos = key if ctg_name is not None else key[1]
        contig = ctg_name if ctg_name is not None else key[0]
        pass_bed_region = len(fp_bed_tree) == 0 or is_region_in(tree=fp_bed_tree,
                                                              contig_name=contig,
                                                              region_start=pos - 1,
                                                              region_end=pos)
        if not pass_bed_region or pos in truth_set:
            pass
            continue

        truth_ref_base = vcf_infos.reference_bases
        truth_alt_base = vcf_infos.alternate_bases[0]
        truth_genotype = vcf_infos.genotype
        if truth_genotype == (0, 0):
            continue
        is_snp_truth = len(truth_ref_base) == 1 and len(truth_alt_base) == 1
        is_ins_truth = len(truth_ref_base) < len(truth_alt_base)
        is_del_truth = len(truth_ref_base) > len(truth_alt_base)

        fn_snp = fn_snp + 1 if is_snp_truth else fn_snp
        fn_ins = fn_ins + 1 if is_ins_truth else fn_ins
        fn_del = fn_del + 1 if is_del_truth else fn_del

        if fn_snp:
            fn_set.add(pos)
    pos_intersection = len(set(truth_variant_dict.keys()).intersection(set(input_variant_dict.keys())))
    print (pos_intersection, len(fp_set), len(fn_set), len(fp_fn_set), len(tp_set), len(fp_set.intersection(fn_set)))
    if output_dir is not None:
        if not os.path.exists(output_dir):
            subprocess.run("mkdir -p {}".format(output_dir), shell=True)
        candidate_types = ['fp', 'fn', 'fp_fn', 'tp']
        variant_sets = [fp_set, fn_set, fp_fn_set, tp_set]
        for vcf_type, variant_set in zip(candidate_types, variant_sets):
            vcf_fn = os.path.join(output_dir, '{}.vcf'.format(vcf_type))
            vcf_writer = VcfWriter(vcf_fn=vcf_fn, ctg_name=ctg_name, write_header=False)
            for pos in variant_set:
                if pos in input_variant_dict:
                    vcf_infos = input_variant_dict[pos]
                elif pos in truth_variant_dict:
                    vcf_infos = truth_variant_dict[pos]
                else:
                    continue
                # ref_base = vcf_infos.reference_bases
                # alt_base = vcf_infos.alternate_bases[0]
                # genotype = vcf_infos.genotype_str
                # qual = float(vcf_infos.qual)
                vcf_writer.write_row(row_str=vcf_infos.row_str)
            vcf_writer.close()

    truth_indel = truth_ins + truth_del
    query_indel = query_ins + query_del
    tp_indel = tp_ins + tp_del
    fp_indel = fp_ins + fp_del
    fn_indel = fn_ins + fn_del
    truth_all = truth_snp + truth_indel
    query_all = query_snp + query_indel
    tp_all = tp_snp + tp_indel
    fp_all = fp_snp + fp_indel
    fn_all = fn_snp + fn_indel

    all_pre, all_rec, all_f1 = cal_metrics(tp=tp_all, fp=fp_all, fn=fn_all)
    snp_pre, snp_rec, snp_f1 = cal_metrics(tp=tp_snp, fp=fp_snp, fn=fn_snp)
    indel_pre, indel_rec, indel_f1 = cal_metrics(tp=tp_indel, fp=fp_indel, fn=fn_indel)
    ins_pre, ins_rec, ins_f1 = cal_metrics(tp=tp_ins, fp=fp_ins, fn=fn_ins)
    del_pre, del_rec, del_f1 = cal_metrics(tp=tp_del, fp=fp_del, fn=fn_del)

    # print (tp_snp, tp_ins, tp_del, fp_snp, fp_ins, fp_del, fn_snp, fn_ins, fn_del, fp_snp_truth, fp_ins_truth, fp_del_truth)
    print ((ctg_name + '-' if ctg_name is not None else "") + input_vcf_fn.split('/')[-1])
    print (len(input_variant_dict), len(truth_variant_dict), pos_out_of_bed)

    print (''.join([item.ljust(15) for item in ["type", 'total.truth', 'total.query', 'tp','fp', 'fn', 'precision', 'recall', "f1-score"]]), file=output_file)
    print (''.join([str(item).ljust(15) for item in ["Overall", truth_all, query_all, tp_all, fp_all, fn_all, all_pre, all_rec, all_f1]]), file=output_file)
    print (''.join([str(item).ljust(15) for item in ["SNP", truth_snp, query_snp, tp_snp, fp_snp, fn_snp, snp_pre, snp_rec, snp_f1]]),file=output_file)
    print (''.join([str(item).ljust(15) for item in ["INDEL", truth_indel, query_indel, tp_indel, fp_indel, fn_indel, indel_pre, indel_rec, indel_f1]]), file=output_file)
    print (''.join([str(item).ljust(15) for item in ["INS", truth_ins, query_ins, tp_ins, fp_ins, fn_ins, ins_pre, ins_rec, ins_f1]]), file=output_file)
    print (''.join([str(item).ljust(15) for item in ["DEL", query_del, query_del, tp_del, fp_del, fn_del, del_pre, del_rec, del_f1]]), file=output_file)
    # print('\n', file=output_file)

    if output_fn:
        output_file.close()



def main():
    parser = ArgumentParser(description="Sort a VCF file according to contig name and starting position")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Input directory")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Input vcf filename suffix")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--contigs_fn', type=str, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--skip_genotyping', action='store_true',
                        help='Output VCF filename, required')

    parser.add_argument('--input_filter_tag', type=str, default=None,
                        help='Output VCF filename, required')

    parser.add_argument('--truth_filter_tag', type=str, default=None,
                        help='Output VCF filename, required')

    args = parser.parse_args()
    compare_vcf(args)

if __name__ == "__main__":
    main()
