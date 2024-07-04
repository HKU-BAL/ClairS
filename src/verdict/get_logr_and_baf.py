from argparse import ArgumentParser

import math

import random
from time import time

seed = int(time())
random.seed(seed)


def getBAFsAndLogRs(tumor_allele_counts_file_prefix, normal_allele_counts_file_prefix, alleles_file_prefix, tumor_logr_output_file, tumor_baf_output_file, normal_baf_output_file, sample_name, normal_sample_name):
    allele_map_dict = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
    allele_dict = dict()
    totalTumor_dict = dict()
    tumorBAF_dict = dict()
    totalNormal_dict = dict()
    normalBAF_dict = dict()
    major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X"]]
    for chr in major_contigs_order:
        alleles_file_path = alleles_file_prefix + str(chr) + '.txt'
        alleles = open(alleles_file_path, 'r')
        for idx, allele in enumerate(alleles.readlines()):
            if idx == 0:
                continue
            allele_info = allele.strip().split('\t')
            allele_pos = allele_info[0]
            allele_ref = allele_map_dict[allele_info[1]]
            allele_alt = allele_map_dict[allele_info[2]]
            allele_dict[(chr, allele_pos)] = allele_ref + '\t' + allele_alt

        tumor_allele_counts_file_path = tumor_allele_counts_file_prefix + str(chr) + '.txt'
        tumor_allele_counts = open(tumor_allele_counts_file_path, 'r')
        for idx, allele_count in enumerate(tumor_allele_counts.readlines()):
            if idx == 0:
                continue
            allele_count_info = allele_count.strip().split('\t')
            allele_chr = allele_count_info[0]
            allele_pos = allele_count_info[1]
            if (chr, allele_pos) not in allele_dict.keys():
                continue
            allele_a_count = allele_count_info[2]
            allele_c_count = allele_count_info[3]
            allele_g_count = allele_count_info[4]
            allele_t_count = allele_count_info[5]
            tumor_allele_ref = allele_dict[(chr, allele_pos)].split('\t')[0]
            tumor_allele_alt = allele_dict[(chr, allele_pos)].split('\t')[1]
            tumor_allele_ref_count = 0
            tumor_allele_alt_count = 0
            if tumor_allele_ref == 'A':
                tumor_allele_ref_count = allele_a_count
            elif tumor_allele_ref == 'C':
                tumor_allele_ref_count = allele_c_count
            elif tumor_allele_ref == 'G':
                tumor_allele_ref_count = allele_g_count
            elif tumor_allele_ref == 'T':
                tumor_allele_ref_count = allele_t_count
            if tumor_allele_alt == 'A':
                tumor_allele_alt_count = allele_a_count
            elif tumor_allele_alt == 'C':
                tumor_allele_alt_count = allele_c_count
            elif tumor_allele_alt == 'G':
                tumor_allele_alt_count = allele_g_count
            elif tumor_allele_alt == 'T':
                tumor_allele_alt_count = allele_t_count
            mutCount1 = int(tumor_allele_ref_count)
            mutCount2 = int(tumor_allele_alt_count)
            totalTumor = mutCount1 + mutCount2
            if totalTumor == 0:
                continue
            key = (allele_chr, allele_pos)
            totalTumor_dict[key] = totalTumor
            tumorBAF = random.choice([(mutCount1 / totalTumor), (mutCount2 / totalTumor)])
            tumorBAF_dict[key] = tumorBAF

        if normal_allele_counts_file_prefix is not None:
            normal_allele_counts_file_path = normal_allele_counts_file_prefix + str(chr) + '.txt'
            normal_allele_counts = open(normal_allele_counts_file_path, 'r')
            for idx, allele_count in enumerate(normal_allele_counts.readlines()):
                if idx == 0:
                    continue
                allele_count_info = allele_count.strip().split('\t')
                allele_chr = allele_count_info[0]
                allele_pos = allele_count_info[1]
                if (chr, allele_pos) not in allele_dict.keys():
                    continue
                allele_a_count = allele_count_info[2]
                allele_c_count = allele_count_info[3]
                allele_g_count = allele_count_info[4]
                allele_t_count = allele_count_info[5]
                normal_allele_ref = allele_dict[(chr, allele_pos)].split('\t')[0]
                normal_allele_alt = allele_dict[(chr, allele_pos)].split('\t')[1]
                normal_allele_ref_count = 0
                normal_allele_alt_count = 0
                if normal_allele_ref == 'A':
                    normal_allele_ref_count = allele_a_count
                elif normal_allele_ref == 'C':
                    normal_allele_ref_count = allele_c_count
                elif normal_allele_ref == 'G':
                    normal_allele_ref_count = allele_g_count
                elif normal_allele_ref == 'T':
                    normal_allele_ref_count = allele_t_count
                if normal_allele_alt == 'A':
                    normal_allele_alt_count = allele_a_count
                elif normal_allele_alt == 'C':
                    normal_allele_alt_count = allele_c_count
                elif normal_allele_alt == 'G':
                    normal_allele_alt_count = allele_g_count
                elif normal_allele_alt == 'T':
                    normal_allele_alt_count = allele_t_count
                mutCount1 = int(normal_allele_ref_count)
                mutCount2 = int(normal_allele_alt_count)
                totalNormal = mutCount1 + mutCount2
                if totalNormal < 10:
                    continue
                key = (allele_chr, allele_pos)
                totalNormal_dict[key] = totalNormal
                normalBAF = random.choice([(mutCount1 / totalNormal), (mutCount2 / totalNormal)])
                normalBAF_dict[key] = normalBAF

    if normal_allele_counts_file_prefix is not None:
        common_keys = [k for k in totalTumor_dict if k in totalNormal_dict]

        tumorBAF_dict = {k: tumorBAF_dict[k] for k in common_keys}
        normalBAF_dict = {k: normalBAF_dict[k] for k in common_keys}

        totalTumor_dict_overlap = {k: totalTumor_dict[k] for k in common_keys}
        totalNormal_dict_overlap = {k: totalNormal_dict[k] for k in common_keys}

        totalTumor_values_overlap = totalTumor_dict_overlap.values()
        mean_totalTumor_overlap = sum(totalTumor_values_overlap) / len(totalTumor_values_overlap)

        totalNormal_values_overlap = totalNormal_dict_overlap.values()
        mean_totalNormal_overlap = sum(totalNormal_values_overlap) / len(totalNormal_values_overlap)

        tumorLogR_dict = {key: math.log2((value / mean_totalTumor_overlap) / (totalNormal_dict_overlap[key] / mean_totalNormal_overlap)) for key, value in totalTumor_dict_overlap.items()}
    else:
        totalTumor_values = totalTumor_dict.values()
        mean_totalTumor = sum(totalTumor_values) / len(totalTumor_values)
        tumorLogR_dict = {key: math.log2(value / mean_totalTumor) for key, value in totalTumor_dict.items()}

    output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
    tumor_logr_output = open(tumor_logr_output_file, 'w')
    tumor_logr_output.write(output_header)
    for key, value in tumorLogR_dict.items():
        tumor_logr_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
        tumor_logr_output.write(tumor_logr_string)
    tumor_baf_output = open(tumor_baf_output_file, 'w')
    tumor_baf_output.write(output_header)
    for key, value in tumorBAF_dict.items():
        tumor_baf_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
        tumor_baf_output.write(tumor_baf_string)
    if normal_allele_counts_file_prefix is not None:
        output_header_normal = 'Chromosome' + '\t' + 'Position' + '\t' + normal_sample_name + '\n'
        normal_baf_output = open(normal_baf_output_file, 'w')
        normal_baf_output.write(output_header_normal)
        for key, value in normalBAF_dict.items():
            normal_baf_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
            normal_baf_output.write(normal_baf_string)


def main():
    parser = ArgumentParser(description="Get Sample LogR and BAF")

    parser.add_argument('--tumor_allele_counts_file_prefix', type=str,
                        default=None,
                        help="Prefix of tumor allele count file")

    parser.add_argument('--normal_allele_counts_file_prefix', type=str,
                        default=None,
                        help="Prefix of normal allele count file")

    parser.add_argument('--alleles_file_prefix', type=str,
                        default=None,
                        help="Prefix of 1kG alleles file")

    parser.add_argument('--tumor_logr_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample LogR")

    parser.add_argument('--tumor_baf_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample BAF")

    parser.add_argument('--normal_baf_output_file', type=str,
                        default=None,
                        help="Output path of normal sample BAF")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")

    parser.add_argument('--normal_sample_name', type=str,
                        default="NORMAL_SAMPLE",
                        help="Normal sample name")

    global args
    args = parser.parse_args()

    getBAFsAndLogRs(args.tumor_allele_counts_file_prefix, args.normal_allele_counts_file_prefix, args.alleles_file_prefix, args.tumor_logr_output_file, args.tumor_baf_output_file, args.normal_baf_output_file, args.sample_name, args.normal_sample_name)


if __name__ == "__main__":
    main()