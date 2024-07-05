from argparse import ArgumentParser

import numpy as np


def make_segments(r, b):
    m = np.column_stack((r, b))

    pcf_segments = []
    previousb = -1
    previousr = 1E10
    count = 0

    for i in range(m.shape[0]):
        if m[i, 1] != previousb or m[i, 0] != previousr:
            if count > 0:
                pcf_segments[-1][-1] = count
            count = 1
            pcf_segments.append([m[i, 0], m[i, 1], count])
        else:
            count += 1
        previousb = m[i, 1]
        previousr = m[i, 0]

    if count > 0:
        pcf_segments[-1][-1] = count

    return np.array(pcf_segments, dtype=float)


def create_distance_matrix(segments, gamma, min_ploidy=None, max_ploidy=None, min_purity=None, max_purity=None):
    s = segments

    if min_ploidy is None or max_ploidy is None:
        psi_pos = np.arange(1, 6.05, 0.05)
    else:
        psi_pos = np.arange(min_ploidy - 0.5, max_ploidy + 0.5, 0.05)

    if min_purity is None or max_purity is None:
        rho_pos = np.arange(0.1, 1.06, 0.01)
    else:
        rho_pos = np.arange(round(min_purity, 2), round(max_purity, 2), 0.01)

    d = np.zeros((len(psi_pos), len(rho_pos)))
    dmin = 1E20

    for i, psi in enumerate(psi_pos):
        for j, rho in enumerate(rho_pos):
            nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

            if np.nansum(nA) < np.nansum(nB):
                nMinor = nA
            else:
                nMinor = nB

            d[i, j] = np.nansum(
                np.abs(nMinor - np.maximum(np.round(nMinor), 0)) ** 2 * s[:, 2] * np.where(s[:, 1] == 0.5, 0.05, 1))

    return d


def rle(x):
    n = len(x)
    y = np.array(x[1:] != x[:-1])
    i = np.append(np.where(y), n - 1)
    lengths = np.diff(np.append(-1, i))
    values = x[i]
    return {'lengths': lengths, 'values': values}


def run_ascat(tumor_logr_file, tumor_baf_file, germline_genotypes_file, tumor_logr_segmented_file, tumor_baf_segmented_file, tumor_purity_ploidy_output_file, tumor_cna_output_file, gamma, min_ploidy, max_ploidy, min_purity,
                    max_purity, sample_name):
    tumor_logr = open(tumor_logr_file, 'r')
    tumor_baf = open(tumor_baf_file, 'r')
    germline_genotypes = open(germline_genotypes_file, 'r')
    tumor_logr_segmented = open(tumor_logr_segmented_file, 'r')
    tumor_baf_segmented = open(tumor_baf_segmented_file, 'r')
    tumor_logr_dict = dict()
    tumor_baf_dict = dict()
    germline_genotypes_dict = dict()
    tumor_logr_segmented_dict = dict()
    tumor_baf_segmented_dict = dict()
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
    for idx, germline_genotypes_line in enumerate(germline_genotypes.readlines()):
        if idx == 0:
            continue
        germline_genotypes_info = germline_genotypes_line.strip().split('\t')
        chr = germline_genotypes_info[0]
        pos = germline_genotypes_info[1]
        gg = germline_genotypes_info[2]
        key = (str(chr), str(pos))
        germline_genotypes_dict[key] = gg
    for idx, tumor_logr_segmented_line in enumerate(tumor_logr_segmented.readlines()):
        if idx == 0:
            continue
        tumor_logr_segmented_info = tumor_logr_segmented_line.strip().split('\t')
        chr = tumor_logr_segmented_info[0]
        pos = tumor_logr_segmented_info[1]
        logr = tumor_logr_segmented_info[2]
        key = (str(chr), str(pos))
        tumor_logr_segmented_dict[key] = logr
    for idx, tumor_baf_segmented_line in enumerate(tumor_baf_segmented.readlines()):
        if idx == 0:
            continue
        tumor_baf_segmented_info = tumor_baf_segmented_line.strip().split('\t')
        chr = tumor_baf_segmented_info[0]
        pos = tumor_baf_segmented_info[1]
        baf = tumor_baf_segmented_info[2]
        key = (str(chr), str(pos))
        tumor_baf_segmented_dict[key] = baf

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

    tumor_baf_ori = np.array(list(tumor_baf_dict.values()))

    germline_genotypes_values = np.array(list(germline_genotypes_dict.values()))
    het = germline_genotypes_values == 'False'
    het_indices = np.where(het)[0]

    b = np.array([float(v) for v in tumor_baf_segmented_dict.values()])
    r = np.array([float(v) for v in tumor_logr_segmented_dict.values()])[het_indices]

    r_ori = np.array([float(v) for v in tumor_logr_segmented_dict.values()])

    s = make_segments(r, b)
    d = create_distance_matrix(s, gamma, min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity,
                               max_purity=max_purity)

    TheoretMaxdist = np.sum(0.25 * s[:, 2] * np.where(s[:, 1] == 0.5, 0.05, 1))

    nonaberrant = False
    MINABB = 0.03
    MINABBREGION = 0.005

    percentAbb = np.sum(np.where(s[:, 1] == 0.5, 0, 1) * s[:, 2]) / np.sum(s[:, 2])
    maxsegAbb = np.max(np.where(s[:, 1] == 0.5, 0, s[:, 2])) / np.sum(s[:, 2])
    if percentAbb <= MINABB and maxsegAbb <= MINABBREGION:
        nonaberrant = True

    MINPLOIDY = min_ploidy
    MAXPLOIDY = max_ploidy
    MINRHO = 0.2
    MINGOODNESSOFFIT = 60
    MINPERCZERO = 0.02
    MINPERCZEROABB = 0.1
    MINPERCODDEVEN = 0.05
    MINPLOIDYSTRICT = 1.7
    MAXPLOIDYSTRICT = 2.3

    localmin = []
    optima = []

    psi_values = np.arange(1.05, 6.05, 0.05)

    rho_values = np.arange(0.11, 1.06, 0.01)
    rho_values = np.round(rho_values, 2)

    for i in range(3, d.shape[0] - 3):
        for j in range(3, d.shape[1] - 3):
            m = d[i, j]
            seld = d[i - 3:i + 4, j - 3:j + 4]
            seld[3, 3] = np.max(seld)
            if np.min(seld) > m:
                psi = psi_values[i]
                rho = rho_values[j]

                nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
                nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

                ploidy = np.sum((nA + nB) * s[:, 2]) / np.sum(s[:, 2])

                percentzero = (np.sum((np.round(nA) == 0) * s[:, 2]) + np.sum((np.round(nB) == 0) * s[:, 2])) / np.sum(
                    s[:, 2])

                goodnessOfFit = (1 - m / TheoretMaxdist) * 100

                if (not nonaberrant and MINPLOIDY < ploidy < MAXPLOIDY and
                        rho >= MINRHO and goodnessOfFit > MINGOODNESSOFFIT and
                        percentzero > MINPERCZERO):
                    optima.append([m, i, j, ploidy, goodnessOfFit])
                    localmin.append(m)

    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if len(optima) == 0 and MINPLOIDY < MAXPLOIDYSTRICT and MAXPLOIDY > MINPLOIDYSTRICT:
        for i in range(3, d.shape[0] - 3):
            for j in range(3, d.shape[1] - 3):
                m = d[i, j]
                seld = d[i - 3:i + 4, j - 3:j + 4]
                seld[3, 3] = np.max(seld)
                if np.min(seld) > m:
                    psi = psi_values[i]
                    rho = rho_values[j]
                    nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
                    nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

                    ploidy = np.sum((nA + nB) * s[:, 2]) / np.sum(s[:, 2])

                    perczeroAbb = (np.sum((np.round(nA) == 0) * s[:, 2] * (s[:, 1] != 0.5)) +
                                   np.sum((np.round(nB) == 0) * s[:, 2] * (s[:, 1] != 0.5))) / \
                                  np.sum(s[:, 2] * (s[:, 1] != 0.5))

                    # Handle the case where BAF is a flat line at 0.5
                    if np.isnan(perczeroAbb):
                        perczeroAbb = 0

                    goodnessOfFit = (1 - m / TheoretMaxdist) * 100

                    if (MINPLOIDYSTRICT < ploidy < MAXPLOIDYSTRICT and
                            rho >= MINRHO and goodnessOfFit > MINGOODNESSOFFIT and
                            perczeroAbb > MINPERCZEROABB):
                        optima.append([m, i, j, ploidy, goodnessOfFit])
                        localmin.append(m)

    # if still no solution, allow solutions with 100% aberrant cells (include the borders with rho = 1), but in first instance, keep the percentzero > 0.01 filter
    if len(optima) == 0:
        # Include borders
        cold = np.where(rho_values > 1)[0]
        d[:, cold] = 1E20

        for i in range(3, d.shape[0] - 3):
            for j in range(3, d.shape[1] - 3):
                m = d[i, j]
                seld = d[i - 3:i + 4, j - 3:j + 4]
                seld[3, 3] = np.max(seld)
                if np.min(seld) > m:
                    psi = psi_values[i]
                    rho = rho_values[j]
                    nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
                    nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

                    ploidy = np.sum((nA + nB) * s[:, 2]) / np.sum(s[:, 2])

                    percentzero = (np.sum((np.round(nA) == 0) * s[:, 2]) + np.sum(
                        (np.round(nB) == 0) * s[:, 2])) / np.sum(s[:, 2])
                    percOddEven = np.sum(((np.round(nA) % 2 == 0) & (np.round(nB) % 2 == 1) | (
                                np.round(nA) % 2 == 1) & (np.round(nB) % 2 == 0)) * s[:, 2]) / np.sum(s[:, 2])
                    perczeroAbb = (np.sum((np.round(nA) == 0) * s[:, 2] * (s[:, 1] != 0.5)) +
                                   np.sum((np.round(nB) == 0) * s[:, 2] * (s[:, 1] != 0.5))) / \
                                  np.sum(s[:, 2] * (s[:, 1] != 0.5))

                    if np.isnan(perczeroAbb):
                        perczeroAbb = 0

                    goodnessOfFit = (1 - m / TheoretMaxdist) * 100

                    if (not nonaberrant and MINPLOIDY < ploidy < MAXPLOIDY and
                            rho >= MINRHO and goodnessOfFit > MINGOODNESSOFFIT and
                            (
                                    perczeroAbb > MINPERCZEROABB or percentzero > MINPERCZERO or percOddEven > MINPERCODDEVEN)):
                        optima.append([m, i, j, ploidy, goodnessOfFit])
                        localmin.append(m)

    # if still no solution, drop the percentzero > MINPERCENTZERO filter, but strict ploidy borders
    if len(optima) == 0 and MINPLOIDY < MAXPLOIDYSTRICT and MAXPLOIDY > MINPLOIDYSTRICT:
        for i in range(3, d.shape[0] - 3):
            for j in range(3, d.shape[1] - 3):
                m = d[i, j]
                seld = d[i - 3:i + 4, j - 3:j + 4]
                seld[3, 3] = np.max(seld)
                if np.min(seld) > m:
                    psi = psi_values[i]
                    rho = rho_values[j]
                    nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
                    nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

                    ploidy = np.sum((nA + nB) * s[:, 2]) / np.sum(s[:, 2])

                    perczeroAbb = (np.sum((np.round(nA) == 0) * s[:, 2] * (s[:, 1] != 0.5)) +
                                   np.sum((np.round(nB) == 0) * s[:, 2] * (s[:, 1] != 0.5))) / \
                                  np.sum(s[:, 2] * (s[:, 1] != 0.5))

                    # Handle the case where BAF is a flat line at 0.5
                    if np.isnan(perczeroAbb):
                        perczeroAbb = 0

                    goodnessOfFit = (1 - m / TheoretMaxdist) * 100

                    if (MINPLOIDYSTRICT < ploidy < MAXPLOIDYSTRICT and
                            rho >= MINRHO and goodnessOfFit > MINGOODNESSOFFIT):
                        optima.append([m, i, j, ploidy, goodnessOfFit])
                        localmin.append(m)

    if optima:
        optlim = np.min(localmin)
        for opt in optima:
            if opt[0] == optlim:
                psi_opt1 = psi_values[int(opt[1])]
                rho_opt1 = rho_values[int(opt[2])]
                if rho_opt1 > 1:
                    rho_opt1 = 1
                ploidy_opt1 = opt[3]
                goodnessOfFit_opt1 = opt[4]
                break

    if optima:
        rho = rho_opt1
        psi = psi_opt1

        SNPposhet = het_indices

        diploidprobes = np.full(len(SNPposhet), True, dtype=bool)

        nAfull = np.where(diploidprobes,
                          (rho - 1 - (b - 1) * 2 ** (r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
                          (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2 ** (r / gamma)) / rho)

        nBfull = np.where(diploidprobes,
                          (rho - 1 + b * 2 ** (r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
                          0)

        nA = np.maximum(np.round(nAfull), 0)
        nB = np.maximum(np.round(nBfull), 0)

        diploidprobes = np.full(len(germline_genotypes_values), True, dtype=bool)

        tlr2 = rle(r_ori)

        tlrstart = np.cumsum(np.concatenate(([0], tlr2['lengths'])))[:-1]
        tlrend = np.cumsum(tlr2['lengths']) - 1
        tlr = tlr2['values']

        seg = []
        seg_raw = []
        for i in range(len(tlr)):
            logR = tlr[i]
            start = tlrstart[i]
            end = tlrend[i]

            baf_slice = np.where((het_indices > start) & (het_indices < end + 1))[0]
            if len(baf_slice) > 0:
                bafke = b[baf_slice][0]
            else:
                baf_slice = np.where((het_indices > start - 100) & (het_indices < end + 1 + 100))[0]
                bafke = b[baf_slice][0]

            nAraw = np.where(diploidprobes[start],
                             (rho - 1 - (bafke - 1) * 2 ** (logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
                             (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2 ** (logR / gamma)) / rho)
            nBraw = np.where(diploidprobes[start],
                             (rho - 1 + bafke * 2 ** (logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho, 0)

            # Correct for negative values
            if nAraw + nBraw < 0:
                nAraw, nBraw = 0, 0
            elif nAraw < 0:
                nBraw += nAraw
                nAraw = 0
            elif nBraw < 0:
                nAraw += nBraw
                nBraw = 0

            # Handle odd copy numbers
            limitround = 0.5
            nA = np.where(bafke == 0.5,
                          np.where(nAraw + nBraw > np.round(nAraw) + np.round(nBraw) + limitround,
                                   np.round(nAraw) + 1,
                                   np.where(nAraw + nBraw < np.round(nAraw) + np.round(nBraw) - limitround,
                                            np.round(nAraw),
                                            np.round(nAraw))),
                          np.round(nAraw))
            nB = np.where(bafke == 0.5,
                          np.where(nAraw + nBraw > np.round(nAraw) + np.round(nBraw) + limitround,
                                   np.round(nBraw),
                                   np.where(nAraw + nBraw < np.round(nAraw) + np.round(nBraw) - limitround,
                                            np.round(nBraw) - 1,
                                            np.round(nBraw))),
                          np.round(nBraw))

            seg.append([start, end, int(nA), int(nB)])
            seg_raw.append([start, end, int(nA), int(nB), float(nAraw), float(nBraw)])

        seg = np.array(seg)
        seg_raw = np.array(seg_raw)

        for _ in range(20):
            seg2 = seg.copy()
            new_seg = []
            skipnext = False
            for i in range(len(seg2)):
                if not skipnext:
                    if (i != len(seg2) - 1 and seg2[i,2] == seg2[i+1,2] and seg2[i,3] == seg2[i+1,3]):
                        segline = [seg2[i,0], seg2[i+1,1], seg2[i,2], seg2[i,3]]
                        skipnext = True
                    else:
                        segline = seg2[i]
                    new_seg.append(segline)
                else:
                    skipnext = False
            seg = np.array(new_seg)

        nMajor = np.zeros(len(r_ori))
        nMinor = np.zeros(len(r_ori))

        for row in seg:
            start, end, nA, nB = row
            nMajor[int(start):int(end)+1] = nA
            nMinor[int(start):int(end)+1] = nB

        n1all = np.zeros(len(r_ori))
        n2all = np.zeros(len(r_ori))

        homo = germline_genotypes_values == 'True'
        homo_indices = np.where(homo)[0]

        heteroprobes_indices = het_indices

        n1all[heteroprobes_indices] = np.where(tumor_baf_ori[heteroprobes_indices].astype(float) <= 0.5, nMajor[heteroprobes_indices], nMinor[heteroprobes_indices])
        n2all[heteroprobes_indices] = np.where(tumor_baf_ori[heteroprobes_indices].astype(float) > 0.5, nMajor[heteroprobes_indices], nMinor[heteroprobes_indices])

        homoprobes_indices = homo_indices

        n1all[homoprobes_indices] = np.where(tumor_baf_ori[homoprobes_indices].astype(float) <= 0.5, nMajor[homoprobes_indices] + nMinor[homoprobes_indices], 0)
        n2all[homoprobes_indices] = np.where(tumor_baf_ori[homoprobes_indices].astype(float) > 0.5, nMajor[homoprobes_indices] + nMinor[homoprobes_indices], 0)

        rho = rho_opt1
        psi = psi_opt1
        goodnessOfFit = goodnessOfFit_opt1
        nonaberrant = nonaberrant
        nA = n1all
        nB = n2all
        seg = seg
        seg_new = []
        for idx, seg_line in enumerate(seg):
            start_idx = int(seg_line[0]) if idx == 0 else int(seg_line[0]) + 1
            end_idx = int(seg_line[1])
            tumor_baf_dict_keys = list(tumor_baf_dict.keys())
            start_key = tumor_baf_dict_keys[start_idx]
            end_key = tumor_baf_dict_keys[end_idx]
            start_chr = start_key[0]
            end_chr = end_key[0]
            start_pos = start_key[1]
            end_pos = end_key[1]
            seg_line_new = [start_chr, start_pos, end_pos, str(seg_line[2]), str(seg_line[3])]
            seg_new.append(seg_line_new)
        seg_raw = seg_raw
        seg_raw_new = []
        for idx, seg_raw_line in enumerate(seg_raw):
            start_idx = int(seg_raw_line[0]) if idx == 0 else int(seg_raw_line[0]) + 1
            end_idx = int(seg_raw_line[1])
            tumor_baf_dict_keys = list(tumor_baf_dict.keys())
            start_key = tumor_baf_dict_keys[start_idx]
            end_key = tumor_baf_dict_keys[end_idx]
            start_chr = start_key[0]
            end_chr = end_key[0]
            start_pos = start_key[1]
            end_pos = end_key[1]
            seg_raw_line_new = [start_chr, start_pos, end_pos, str(seg_raw_line[2]), str(seg_raw_line[3])]
            seg_raw_new.append(seg_raw_line_new)
        distance_matrix = d
        ploidy = np.mean(n1all + n2all)

    else:
        print("Could not find an optimal purity and ploidy value for {}!".format(sample_name))
        rho = None
        psi = None
        goodnessOfFit = None
        nonaberrant = None
        nA = None
        nB = None
        seg = None
        seg_new = None
        seg_raw = None
        seg_raw_new = None
        distance_matrix = None
        ploidy = None

    if rho != None:
        output_header_purity = 'Sample' + '\t' + 'Purity' + '\t' + 'Ploidy' + '\t' + 'GoodnessOfFit' + '\n'
        tumor_purity_ploidy_output = open(tumor_purity_ploidy_output_file, 'w')
        tumor_purity_ploidy_output.write(output_header_purity)
        tumor_purity_ploidy_string = sample_name + '\t' + str(rho) + '\t' + str(ploidy) + '\t' + str(goodnessOfFit) + '\n'
        tumor_purity_ploidy_output.write(tumor_purity_ploidy_string)

        output_header_cna = 'Sample' + '\t' + 'Chromosome' + '\t' + 'StartPosition' + '\t' + 'EndPosition' + '\t' + 'nMajor' + '\t' + 'nMinor' + '\n'
        tumor_cna_output = open(tumor_cna_output_file, 'w')
        tumor_cna_output.write(output_header_cna)
        for seg_line in seg_new:
            tumor_cna_string = sample_name + '\t' + '\t'.join(list(seg_line)) + '\n'
            tumor_cna_output.write(tumor_cna_string)


def main():
    parser = ArgumentParser(description="Run ASCAT")

    parser.add_argument('--tumor_logr_file', type=str,
                        default=None,
                        help="Path of tumor sample LogR")

    parser.add_argument('--tumor_baf_file', type=str,
                        default=None,
                        help="Path of tumor sample BAF")

    parser.add_argument('--germline_genotypes_file', type=str,
                        default=None,
                        help="Path of germline genotypes")

    parser.add_argument('--tumor_logr_segmented_file', type=str,
                        default=None,
                        help="Path of tumor sample PCFed LogR")

    parser.add_argument('--tumor_baf_segmented_file', type=str,
                        default=None,
                        help="Path of tumor sample PCFed BAF")

    parser.add_argument('--tumor_purity_ploidy_output_file', type=str,
                        default=None,
                        help="Output path of estimated tumor sample purity and ploidy")

    parser.add_argument('--tumor_cna_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample CNA file")

    parser.add_argument('--gamma', type=float,
                        default=1.0,
                        help="Value of gamma parameter")

    parser.add_argument('--min_ploidy', type=float,
                        default=1.5,
                        help="Value of min ploidy")

    parser.add_argument('--max_ploidy', type=float,
                        default=5.5,
                        help="Value of max ploidy")

    parser.add_argument('--min_purity', type=float,
                        default=0.1,
                        help="Value of min purity")

    parser.add_argument('--max_purity', type=float,
                        default=1.05,
                        help="Value of max purity")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")

    global args
    args = parser.parse_args()

    run_ascat(args.tumor_logr_file, args.tumor_baf_file, args.germline_genotypes_file, args.tumor_logr_segmented_file, args.tumor_baf_segmented_file, args.tumor_purity_ploidy_output_file, args.tumor_cna_output_file, args.gamma, args.min_ploidy,
              args.max_ploidy, args.min_purity, args.max_purity, args.sample_name)


if __name__ == "__main__":
    main()