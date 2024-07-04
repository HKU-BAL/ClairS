from argparse import ArgumentParser

import numpy as np
import math


def predictGermlineGenotypes(tumor_logr_file, tumor_baf_file, normal_baf_file, germline_genotypes_output_file, maxHomozygous, proportionHetero, proportionHomo, proportionOpen, segmentLength, sample_name):
    if normal_baf_file is None:
        tumor_logr = open(tumor_logr_file, 'r')
        tumor_baf = open(tumor_baf_file, 'r')
        tumor_logr_dict = dict()
        tumor_baf_dict = dict()
        for idx, tumor_logr_line in enumerate(tumor_logr.readlines()):
            if idx == 0:
                continue
            tumor_logr_info = tumor_logr_line.strip().split('\t')
            chr = tumor_logr_info[0]
            pos = tumor_logr_info[1]
            logr = tumor_logr_info[2]
            key = (str(chr), str(pos))
            tumor_logr_dict[key] = logr
        for idx, tumor_baf_line in enumerate(tumor_baf.readlines()):
            if idx == 0:
                continue
            tumor_baf_info = tumor_baf_line.strip().split('\t')
            chr = tumor_baf_info[0]
            pos = tumor_baf_info[1]
            baf = tumor_baf_info[2]
            key = (str(chr), str(pos))
            tumor_baf_dict[key] = baf

        result = []
        current_chr = None
        temp_list = []
        for index, (key, value) in enumerate(tumor_baf_dict.items()):
            chr_key = key[0]

            if chr_key != current_chr:
                if temp_list:
                    result.append(temp_list)
                temp_list = []
                current_chr = chr_key

            temp_list.append(index)

        if temp_list:
            result.append(temp_list)

        tbsam = np.array(list(tumor_baf_dict.values())).astype(float)
        bsm = np.where(tbsam < 0.5, tbsam, 1 - tbsam)

        sorted_bsm = np.sort(bsm)
        index = round(len(bsm) * proportionHomo)
        value = sorted_bsm[index]
        homoLimit = max(value, maxHomozygous)

        Hom = np.where(bsm < homoLimit, True, np.nan)
        Undecided = np.sum(np.isnan(Hom))
        extraHetero = round(min(proportionHetero * len(tbsam), Undecided - proportionOpen * len(tbsam)))

        if extraHetero > 0:
            allProbes = np.arange(len(tbsam))
            nonHomoProbes = allProbes[np.isnan(Hom) | (Hom == False)]

            lowestDist = []

            bsmHNA = bsm.copy()
            bsmHNA[~np.isnan(Hom) & (Hom == True)] = np.nan

            for chr_result in result:
                chrNonHomoProbes = sorted(list(set(nonHomoProbes).intersection(set(chr_result))))

                if len(chrNonHomoProbes) > 5:
                    segmentLength2 = min(len(chrNonHomoProbes) - 1, segmentLength)

                    # Window on the left
                    chrNonHomoProbesStartWindowLeft = np.concatenate((np.repeat(np.nan, segmentLength2), chrNonHomoProbes[:len(chrNonHomoProbes) - segmentLength2]))
                    chrNonHomoProbesEndWindowLeft = np.concatenate(([np.nan], chrNonHomoProbes[:len(chrNonHomoProbes) - 1]))

                    # Window on the right
                    chrNonHomoProbesStartWindowRight = np.concatenate((chrNonHomoProbes[1:], [np.nan]))
                    chrNonHomoProbesEndWindowRight = np.concatenate((chrNonHomoProbes[segmentLength2:], np.repeat(np.nan, segmentLength2)))

                    # Window in the middle
                    middle_segment_length = segmentLength2 // 2
                    chrNonHomoProbesStartWindowMiddle = np.concatenate((np.repeat(np.nan, middle_segment_length), chrNonHomoProbes[:len(chrNonHomoProbes) - middle_segment_length]))
                    chrNonHomoProbesEndWindowMiddle = np.concatenate((chrNonHomoProbes[middle_segment_length:], np.repeat(np.nan, middle_segment_length)))

                    chrLowestDist = []

                    for probeNr in range(len(chrNonHomoProbes)):
                        probe = chrNonHomoProbes[probeNr]

                        start_left = chrNonHomoProbesStartWindowLeft[probeNr]
                        end_left = chrNonHomoProbesEndWindowLeft[probeNr]

                        if not math.isnan(start_left) and not math.isnan(end_left):
                            window_values_left = bsmHNA[int(start_left):int(end_left) + 1]
                            window_values_left = [value for value in window_values_left if not math.isnan(value)]
                            if window_values_left:
                                medianLeft = np.median(window_values_left)
                            else:
                                medianLeft = np.nan
                        else:
                            medianLeft = np.nan

                        start_right = chrNonHomoProbesStartWindowRight[probeNr]
                        end_right = chrNonHomoProbesEndWindowRight[probeNr]

                        if not math.isnan(start_right) and not math.isnan(end_right):
                            window_values_right = bsmHNA[int(start_right):int(end_right) + 1]
                            window_values_right = [value for value in window_values_right if not math.isnan(value)]

                            if window_values_right:
                                medianRight = np.median(window_values_right)
                            else:
                                medianRight = np.nan
                        else:
                            medianRight = np.nan

                        start_middle = chrNonHomoProbesStartWindowMiddle[probeNr]
                        end_middle = chrNonHomoProbesEndWindowMiddle[probeNr]

                        if not math.isnan(start_middle) and not math.isnan(end_middle):
                            window_values_middle_left = bsmHNA[int(start_middle):int(end_left) + 1]
                            window_values_middle_right = bsmHNA[int(start_right):int(end_middle) + 1]

                            window_values_middle_left = [value for value in window_values_middle_left if
                                                         not math.isnan(value)]
                            window_values_middle_right = [value for value in window_values_middle_right if
                                                          not math.isnan(value)]

                            concatenated_values = np.concatenate(
                                (window_values_middle_left, window_values_middle_right))

                            if concatenated_values.size > 0:
                                medianMiddle = np.median(concatenated_values)
                            else:
                                medianMiddle = np.nan
                        else:
                            medianMiddle = np.nan

                        diffs = [np.abs(medianLeft - bsm[probe]),
                                 np.abs(medianRight - bsm[probe]),
                                 np.abs(medianMiddle - bsm[probe])]
                        diffs = [d for d in diffs if not np.isnan(d)]
                        chrLowestDist.insert(probeNr, np.min(diffs) if len(diffs) > 0 else np.inf)
                else:
                    if len(chrNonHomoProbes) > 0:
                        chrLowestDist = [1] * len(chrNonHomoProbes)
                    else:
                        chrLowestDist = []

                lowestDist.extend(chrLowestDist)

            lowestDistUndecided = [lowestDist[i] for i, h in enumerate(Hom[nonHomoProbes]) if np.isnan(h)]
            sorted_indices = np.argsort(lowestDistUndecided)
            sorted_lowestDistUndecided = np.sort(lowestDistUndecided)

            min_length = min(len(sorted_lowestDistUndecided), extraHetero)
            selected_indices = sorted_indices[:min_length]
            selected_indices_ori = [nonHomoProbes[i] for i in selected_indices]

            Hom[selected_indices_ori] = False

        Hom[np.isnan(Hom)] = True
        germline_genotypes_output_dict = {key: Hom[idx] for idx, key in enumerate(tumor_baf_dict.keys())}

        output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
        germline_genotypes_output = open(germline_genotypes_output_file, 'w')
        germline_genotypes_output.write(output_header)
        for key, value in germline_genotypes_output_dict.items():
            if value == 1.0:
                value = "True"
            elif value == 0.0:
                value = "False"
            germline_genotypes_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
            germline_genotypes_output.write(germline_genotypes_string)
    else:
        normal_baf = open(normal_baf_file, 'r')
        normal_baf_dict = dict()
        for idx, normal_baf_line in enumerate(normal_baf.readlines()):
            if idx == 0:
                continue
            normal_baf_info = normal_baf_line.strip().split('\t')
            chr = normal_baf_info[0]
            pos = normal_baf_info[1]
            baf = normal_baf_info[2]
            key = (str(chr), str(pos))
            normal_baf_dict[key] = baf
        normal_baf_values = np.array(list(normal_baf_dict.values()), dtype=float)
        Hom = np.where((normal_baf_values < 0.3) | (normal_baf_values > 0.7), 'True', 'False')
        germline_genotypes_output_dict = {key: Hom[idx] for idx, key in enumerate(normal_baf_dict.keys())}
        output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
        germline_genotypes_output = open(germline_genotypes_output_file, 'w')
        germline_genotypes_output.write(output_header)
        for key, value in germline_genotypes_output_dict.items():
            germline_genotypes_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
            germline_genotypes_output.write(germline_genotypes_string)


def main():
    parser = ArgumentParser(description="Predict Germline Genotypes")

    parser.add_argument('--tumor_logr_file', type=str,
                        default=None,
                        help="Path of tumor sample LogR")

    parser.add_argument('--tumor_baf_file', type=str,
                        default=None,
                        help="Pclear"
                             "Path of tumor sample BAF")

    parser.add_argument('--normal_baf_file', type=str,
                        default=None,
                        help="Output path of normal sample BAF")

    parser.add_argument('--germline_genotypes_output_file', type=str,
                        default=None,
                        help="Output path of germline genotypes")

    parser.add_argument('--maxHomozygous', type=float,
                        default=0.02,
                        help="Value of max homozygous")

    parser.add_argument('--proportionHetero', type=float,
                        default=0.30,
                        help="Proportion of hetero")

    parser.add_argument('--proportionHomo', type=float,
                        default=0.65,
                        help="Proportion of homo")

    parser.add_argument('--proportionOpen', type=float,
                        default=0.03,
                        help="Proportion of open")

    parser.add_argument('--segmentLength', type=int,
                        default=100,
                        help="Segment length")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")

    global args
    args = parser.parse_args()

    predictGermlineGenotypes(args.tumor_logr_file, args.tumor_baf_file, args.normal_baf_file, args.germline_genotypes_output_file, args.maxHomozygous, args.proportionHetero, args.proportionHomo,
                             args.proportionOpen, args.segmentLength, args.sample_name)


if __name__ == "__main__":
    main()