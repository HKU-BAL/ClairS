from argparse import ArgumentParser

import numpy as np
import math
from scipy.ndimage import median_filter

import random
from time import time

seed = int(time())
random.seed(seed)


def predictGermlineHomozygousStretches(chr, hom):
    homsam = hom
    num_hom = np.sum(homsam == 'True')
    total_num = len(homsam)
    perchom = num_hom / total_num
    if perchom == 0.0:
        homthres = 0
    elif perchom == 1.0:
        homthres = 1
    else:
        homthres = math.ceil(math.log(0.001, perchom))
    allhprobes = []
    for chrke, chrom in enumerate(chr):
        hschr = homsam[chrom]
        hprobes = []
        for probe, value in enumerate(hschr):
            if value == 'True':
                hprobes.append(probe)
            else:
                if len(hprobes) >= homthres:
                    allhprobes.append([chrke, chrom[min(hprobes)], chrom[max(hprobes)]])
                hprobes = []
        # if the last probe is homozygous, this is not yet accounted for
        if hschr[-1] == 'True':
            if len(hprobes) >= homthres:
                allhprobes.append([chrke, chrom[min(hprobes)], chrom[max(hprobes)]])

    if not allhprobes:
        allhprobes = [[0, 0, 0]]

    HomoStretches = allhprobes

    return HomoStretches


def fastAspcf(logR, allB, kmin, gamma):
    N = len(logR)
    w = 1000  # w: windowsize
    d = 100

    startw = -d
    stopw = w - d

    nseg = 0
    var2 = 0
    var3 = 0
    breakpts = [0]

    while True:
        from_idx = max(0, startw)
        to_idx = min(stopw, N)
        logRpart = logR[from_idx:to_idx]  # Adjusted for Python's 0-based index
        allBpart = allB[from_idx:to_idx]
        allBflip = allBpart.copy()
        allBflip[allBpart > 0.5] = 1 - allBpart[allBpart > 0.5]

        sd1 = getMad(logRpart)
        sd2 = getMad(allBflip)
        sd3 = getMad(allBpart)

        # Must check that sd1 and sd2 are defined and != 0:
        sd_valid = [not np.isnan(sd1), not np.isnan(sd2), sd1 != 0, sd2 != 0]
        if all(sd_valid):
            # run aspcfpart:
            part_res = aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N,
                                 kmin=kmin, gamma=gamma)
            breakpts_part = part_res['breakpts']
            breakpts_part_np = np.array(breakpts_part)
            last_breakpt = breakpts[-1]
            larger = breakpts_part_np > last_breakpt
            breakpts.extend(breakpts_part_np[larger])
            var2 += sd2 ** 2
            var3 += sd3 ** 2
            nseg += 1

        if stopw < N + d:
            startw = min(stopw - 2 * d + 1, N - 2 * d)
            stopw = startw + w
        else:
            break

    breakpts = list(np.unique(breakpts + [N]))
    if nseg == 0:
        nseg = 1  # just in case the sd-test never passes.
    sd2 = np.sqrt(var2 / nseg)
    sd3 = np.sqrt(var3 / nseg)

    # On each segment calculate mean of unflipped B allele data
    frst = np.array(breakpts[:-1]) + 1
    last = np.array(breakpts[1:])
    nseg = len(frst)

    yhat1 = np.full(N, np.nan)
    yhat2 = np.full(N, np.nan)

    for i in range(nseg):
        yhat1[frst[i] - 1:last[i]] = np.mean(logR[frst[i] - 1:last[i]])
        yi2 = allB[frst[i] - 1:last[i]]
        # Center data around zero (by subtracting 0.5) and estimate mean
        if len(yi2) == 0:
            mu = 0
        else:
            mu = np.mean(np.abs(yi2 - 0.5))

        if np.sqrt(sd2 ** 2 + mu ** 2) < 2 * sd2:
            mu = 0

        yhat2[frst[i] - 1:last[i]] = mu + 0.5

    return {'yhat1': yhat1, 'yhat2': yhat2}


def aspcfpart(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma):
    from_idx = max(0, a)
    usefrom = max(0, a + d)
    useto = min(N, b - d)

    N = len(logRpart)
    y1 = logRpart
    y2 = allBflip

    # Check that vectors are long enough to run algorithm:
    if N < 2 * kmin:
        breakpts = [0]
        return {'breakpts': breakpts}

    initSum1 = np.sum(y1[:kmin])
    initKvad1 = np.sum(y1[:kmin] ** 2)
    initAve1 = initSum1 / kmin
    initSum2 = np.sum(y2[:kmin])
    initKvad2 = np.sum(y2[:kmin] ** 2)
    initAve2 = initSum2 / kmin

    # Define vector of best costs
    bestCost = np.zeros(N)
    cost1 = (initKvad1 - initSum1 * initAve1) / sd1 ** 2
    cost2 = (initKvad2 - initSum2 * initAve2) / sd2 ** 2
    bestCost[kmin - 1] = cost1 + cost2

    # Define vector of best splits
    bestSplit = np.zeros(N, dtype=int)

    # Define vector of best averages
    bestAver1 = np.zeros(N)
    bestAver2 = np.zeros(N)
    bestAver1[kmin - 1] = initAve1
    bestAver2[kmin - 1] = initAve2

    # Initialize
    Sum1 = np.zeros(N)
    Sum2 = np.zeros(N)
    Kvad1 = np.zeros(N)
    Kvad2 = np.zeros(N)
    Aver1 = np.zeros(N)
    Aver2 = np.zeros(N)
    Cost = np.zeros(N)

    kminP1 = kmin + 1
    for k in range(kminP1, 2 * kmin):
        Sum1[kminP1 - 1:k] += y1[k - 1]
        Aver1[kminP1 - 1:k] = Sum1[kminP1 - 1:k] / np.arange(k - kmin, 0, -1)
        Kvad1[kminP1 - 1:k] += y1[k - 1] ** 2
        Sum2[kminP1 - 1:k] += y2[k - 1]
        Aver2[kminP1 - 1:k] = Sum2[kminP1 - 1:k] / np.arange(k - kmin, 0, -1)
        Kvad2[kminP1 - 1:k] += y2[k - 1] ** 2

        bestAver1[k - 1] = (initSum1 + Sum1[kminP1 - 1]) / k
        bestAver2[k - 1] = (initSum2 + Sum2[kminP1 - 1]) / k
        cost1 = ((initKvad1 + Kvad1[kminP1 - 1]) - k * bestAver1[k - 1] ** 2) / sd1 ** 2
        cost2 = ((initKvad2 + Kvad2[kminP1 - 1]) - k * bestAver2[k - 1] ** 2) / sd2 ** 2

        bestCost[k - 1] = cost1 + cost2

    for n in range(2 * kmin, N + 1):
        nMkminP1 = n - kmin + 1

        Sum1[kminP1 - 1:n] += y1[n - 1]
        Aver1[kminP1 - 1:n] = Sum1[kminP1 - 1:n] / np.arange(n - kmin, 0, -1)
        Kvad1[kminP1 - 1:n] += y1[n - 1] ** 2

        cost1 = (Kvad1[kminP1 - 1:nMkminP1] - Sum1[kminP1 - 1:nMkminP1] * Aver1[kminP1 - 1:nMkminP1]) / sd1 ** 2

        Sum2[kminP1 - 1:n] += y2[n - 1]
        Aver2[kminP1 - 1:n] = Sum2[kminP1 - 1:n] / np.arange(n - kmin, 0, -1)
        Kvad2[kminP1 - 1:n] += y2[n - 1] ** 2
        cost2 = (Kvad2[kminP1 - 1:nMkminP1] - Sum2[kminP1 - 1:nMkminP1] * Aver2[kminP1 - 1:nMkminP1]) / sd2 ** 2

        Cost[kminP1 - 1:nMkminP1] = bestCost[kmin - 1:n - kmin] + cost1 + cost2

        Pos = np.argmin(Cost[kminP1 - 1:nMkminP1]) + kmin
        cost = Cost[Pos - 1] + gamma

        aver1 = Aver1[Pos - 1]
        aver2 = Aver2[Pos - 1]
        totAver1 = (Sum1[kminP1 - 1] + initSum1) / n
        totCost1 = ((Kvad1[kminP1 - 1] + initKvad1) - n * totAver1 ** 2) / sd1 ** 2
        totAver2 = (Sum2[kminP1 - 1] + initSum2) / n
        totCost2 = ((Kvad2[kminP1 - 1] + initKvad2) - n * totAver2 ** 2) / sd2 ** 2
        totCost = totCost1 + totCost2

        if totCost < cost:
            Pos = 1
            cost = totCost
            aver1 = totAver1
            aver2 = totAver2

        bestCost[n - 1] = cost
        bestAver1[n - 1] = aver1
        bestAver2[n - 1] = aver2
        bestSplit[n - 1] = Pos - 1

    # Trace back
    n = N
    breakpts = [n]
    while n > 0:
        breakpts.append(bestSplit[n - 1])
        n = bestSplit[n - 1]

    breakpts = np.array(breakpts) + from_idx - 1
    breakpts = breakpts[(breakpts >= usefrom) & (breakpts <= useto)].tolist()

    return {'breakpts': breakpts}


def getMad(x, k=25):
    # Remove observations that are equal to zero
    x = x[x != 0]

    run_median = medianFilter(x, k)

    # Calculate differences
    dif = x - run_median

    # Calculate MAD (Median Absolute Deviation)
    mad = np.median(np.abs(dif - np.median(dif)))

    return mad


def exactPcf(y, kmin, gamma):
    """
    Implementation of exact PCF by Potts-filtering.

    Parameters:
    y (array-like): Input array of (log2) copy numbers.
    kmin (int): Minimal length of plateaus.
    gamma (float): Penalty for each discontinuity.

    Returns:
    np.ndarray: Output array of filtered values.
    """
    N = len(y)
    y = np.array(y)
    yhat = np.zeros(N)

    if N < 2 * kmin:
        yhat[:] = np.mean(y)
        return yhat

    init_sum = np.sum(y[:kmin])
    init_kvad = np.sum(y[:kmin] ** 2)
    init_ave = init_sum / kmin

    best_cost = np.zeros(N)
    best_cost[kmin - 1] = init_kvad - init_sum * init_ave

    best_split = np.zeros(N, dtype=int)
    best_aver = np.zeros(N)
    best_aver[kmin - 1] = init_ave

    Sum = np.zeros(N)
    Kvad = np.zeros(N)
    Aver = np.zeros(N)
    Cost = np.zeros(N)

    kminP1 = kmin + 1

    for k in range(kminP1, 2 * kmin):
        Sum[kminP1 - 1:k] += y[k - 1]
        Aver[kminP1 - 1:k] = Sum[kminP1 - 1:k] / np.arange(k - kmin, 0, -1)
        Kvad[kminP1 - 1:k] += y[k - 1] ** 2
        best_aver[k - 1] = (init_sum + Sum[kminP1 - 1]) / k
        best_cost[k - 1] = (init_kvad + Kvad[kminP1 - 1]) - k * best_aver[k - 1] ** 2

    for n in range(2 * kmin, N + 1):
        yn = y[n - 1]
        yn2 = yn ** 2
        Sum[kminP1 - 1:n] += yn
        Aver[kminP1 - 1:n] = Sum[kminP1 - 1:n] / np.arange(n - kmin, 0, -1)
        Kvad[kminP1 - 1:n] += yn2

        nMkminP1 = n - kmin + 1
        Cost[kminP1 - 1:nMkminP1] = best_cost[kmin - 1:n - kmin] + Kvad[kminP1 - 1:nMkminP1] - Sum[kminP1 - 1:nMkminP1] * Aver[kminP1 - 1:nMkminP1] + gamma

        Pos = np.argmin(Cost[kminP1 - 1:nMkminP1]) + kmin
        cost = Cost[Pos - 1]
        aver = Aver[Pos - 1]

        tot_aver = (Sum[kminP1 - 1] + init_sum) / n
        tot_cost = (Kvad[kminP1 - 1] + init_kvad) - n * tot_aver ** 2

        if tot_cost < cost:
            Pos = 1
            cost = tot_cost
            aver = tot_aver

        best_cost[n - 1] = cost
        best_aver[n - 1] = aver
        best_split[n - 1] = Pos - 1

    n = N
    while n > 0:
        yhat[best_split[n - 1]:n] = best_aver[n - 1]
        n = best_split[n - 1]

    return {'yhat': yhat}


def madWins(x, tau, k):
    xhat = medianFilter(x, k)

    # Calculate differences
    d = x - xhat

    # Calculate MAD (Median Absolute Deviation)
    mad = np.median(np.abs(d - np.median(d)))

    # Calculate threshold z
    z = tau * mad

    # Apply psi function
    xwin = xhat + psi(d, z)

    # Identify outliers
    outliers = np.zeros(len(x))
    outliers[x > xwin] = 1
    outliers[x < xwin] = -1

    return {'ywin': xwin, 'sdev': mad, 'outliers': outliers}


def medianFilter(x, k):
    n = len(x)
    filtWidth = 2 * k + 1

    # Ensure filtWidth does not exceed n
    if filtWidth > n:
        if n == 0:
            filtWidth = 1
        elif n % 2 == 0:
            # Ensure filtWidth is odd
            filtWidth = n - 1
        else:
            filtWidth = n

    # Apply running median filter
    runMedian = median_filter(x, size=filtWidth, mode='reflect')

    return runMedian


def psi(x, z):
    """
    Clip the values in the array x to the range [-z, z].

    Parameters:
    x (array-like): Input array.
    z (float): Threshold value.

    Returns:
    numpy.ndarray: Clipped array.
    """
    xwin = np.copy(x)
    xwin[x < -z] = -z
    xwin[x > z] = z
    return xwin


def fillNA(x, zeroIsNA=False):
    """
    Fill NaN values in the input array.

    Args:
        x (numpy.ndarray): Input array.
        zeroIsNA (bool): If True, treat zeros as NaN.

    Returns:
        numpy.ndarray: Array with NaN values filled.
    """
    x_filled = x.copy()
    if zeroIsNA:
        x_filled[x_filled == 0] = np.nan

    # Fill NaN values using linear interpolation
    idx = np.where(~np.isnan(x_filled))[0]
    x_filled[np.isnan(x_filled)] = np.interp(np.where(np.isnan(x_filled))[0], idx, x_filled[idx])

    return x_filled


def rle_lengths(arr):
    """Run-length encoding to get lengths of consecutive values."""
    n = len(arr)
    if n == 0:
        return np.array([], dtype=int)
    else:
        y = np.array(arr[1:] != arr[:-1])
        i = np.append(np.where(y)[0], n - 1)
        return np.diff(np.append(-1, i))


def aspcf(tumor_logr_file, tumor_baf_file, germline_genotypes_file, tumor_logr_pcfed_output_file, tumor_baf_pcfed_output_file, penalty, sample_name):
    tumor_logr = open(tumor_logr_file, 'r')
    tumor_baf = open(tumor_baf_file, 'r')
    germline_genotypes = open(germline_genotypes_file, 'r')
    tumor_logr_dict = dict()
    tumor_baf_dict = dict()
    germline_genotypes_dict = dict()
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

    germline_genotypes_values = np.array(list(germline_genotypes_dict.values()))

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

    ghs = predictGermlineHomozygousStretches(result, germline_genotypes_values)

    segmentlengths = sorted(list({penalty, 70, 100, 140}))
    segmentlengths = [l for l in segmentlengths if l >= penalty]

    logRPCFed = np.array([])
    bafPCFed = np.array([])

    for segmentlength in segmentlengths:
        logRPCFed = np.array([])
        bafPCFed = np.array([])
        tbsam = np.array(list(tumor_baf_dict.values())).astype(float)
        logr = np.array(list(tumor_logr_dict.values())).astype(float)
        homosam = np.array(list(germline_genotypes_dict.values()))
        for chrke, chrom in enumerate(result):
            lr = logr[chrom]
            lrwins = madWins(lr, 2.5, 25)['ywin']
            baf = tbsam[chrom]
            homo = homosam[chrom]
            Select_het = homo == 'False'
            bafsel = baf[Select_het]
            bafselwinsmirrored = madWins(np.where(bafsel > 0.5, bafsel, 1 - bafsel), 2.5, 25)['ywin']
            bafselwins = np.where(bafsel > 0.5, bafselwinsmirrored, 1 - bafselwinsmirrored)
            select_het_indices = np.where(Select_het)[0]
            logRaveraged = None

            if len(select_het_indices) != 0:
                averageIndices = np.concatenate(([0], (select_het_indices[:-1] + select_het_indices[1:]) / 2, [len(lr)]))
                startindices = np.ceil(averageIndices[:-1]).astype(int)
                endindices = np.floor(averageIndices[1:]).astype(int)

                if len(select_het_indices) == 1:
                    startindices = [0]
                    endindices = [len(lr) - 1]

                logRaveraged = np.full(len(select_het_indices), np.nan)

                for i in range(len(select_het_indices)):
                    if np.isnan(endindices[i]):
                        endindices[i] = startindices[i]
                    logRaveraged[i] = np.nanmean(lrwins[startindices[i]:endindices[i] + 1])

            if len(logRaveraged) > 0:
                if len(logRaveraged) < 6:
                    logRASPCF = np.full(len(logRaveraged), np.mean(logRaveraged))
                    bafASPCF = np.full(len(logRaveraged), np.mean(bafselwinsmirrored))
                else:
                    PCFed = fastAspcf(logRaveraged, bafselwins, 6, segmentlength)
                    logRASPCF = PCFed['yhat1']
                    bafASPCF = PCFed['yhat2']

                logRc = np.array([], dtype=float)

                for probe in range(len(logRASPCF)):
                    if probe == 0:
                        logRc = np.concatenate((logRc, np.full(select_het_indices[probe], logRASPCF[probe])))
                    elif probe == len(logRASPCF) - 1:
                        logRc = np.concatenate((logRc, np.full(len(lr) - select_het_indices[probe], logRASPCF[probe])))
                    else:
                        start = select_het_indices[probe]
                        end = select_het_indices[probe + 1]
                        interval_length = end - start

                        if logRASPCF[probe] == logRASPCF[probe + 1]:
                            logRc = np.concatenate((logRc, np.full(interval_length, logRASPCF[probe])))
                        else:
                            d = np.array([])
                            for bp in range(interval_length):
                                dis = np.sum(np.abs(lr[start:start + bp] - logRASPCF[probe]))
                                dis += np.sum(np.abs(lr[start + bp + 1:end] - logRASPCF[probe + 1]))
                                d = np.append(d, dis)

                            breakpoint = np.argmin(d)
                            logRc = np.concatenate((logRc, np.full(breakpoint, logRASPCF[probe]),
                                                    np.full(interval_length - breakpoint, logRASPCF[probe + 1])))

                last_length = len(lr) - len(logRc)
                logRc = np.concatenate((logRc, np.full(last_length, logRASPCF[-1])))

                seg = rle_lengths(logRc)

                # Initialize variables
                logRd = np.array([], dtype=float)
                startprobe = 0

                for length in seg:
                    endprobe = startprobe + length
                    level = np.nanmean(lr[startprobe:endprobe])  # Compute mean ignoring NaNs
                    logRd = np.concatenate((logRd, np.full(length, level)))
                    startprobe = endprobe

                logRPCFed = np.concatenate((logRPCFed, logRd))
                bafPCFed = np.concatenate((bafPCFed, bafASPCF))

            else:
                level = np.nanmean(lr)
                reps = len(lr)
                logRPCFed = np.concatenate((logRPCFed, np.full(reps, level)))

            homsegs = [sub_list for sub_list in ghs if sub_list[0] == chrke]
            homsegs = np.array(homsegs)
            startchr = np.min(result[chrke])
            endchr = np.max(result[chrke])

            if homsegs is not None and len(homsegs) != 0:
                for i in range(len(homsegs)):
                    startpos = max(homsegs[i, 1], startchr)
                    endpos = min(homsegs[i, 2], endchr)

                    startpos2 = max(homsegs[i, 1] - 100, startchr)
                    endpos2 = min(homsegs[i, 2] + 100, endchr)

                    startpos3 = max(homsegs[i, 1] - 5, startchr)
                    endpos3 = min(homsegs[i, 2] + 5, endchr)

                    towins = logr[startpos2:endpos2 + 1]
                    winsed = madWins(towins[~np.isnan(towins)], 2.5, 25)['ywin']
                    pcfed = np.full(len(towins), np.nan)
                    pcfed[~np.isnan(towins)] = exactPcf(winsed, 6, int(segmentlength / 4))['yhat']
                    pcfed2 = pcfed[startpos3 - startpos2:endpos3 - startpos2 + 1]

                    if len(pcfed2) != len(logRPCFed[startpos3:endpos3 + 1]):
                        pcfed2 = pcfed2[:len(logRPCFed[startpos3:endpos3 + 1])]
                    dif = np.abs(pcfed2 - logRPCFed[startpos3:endpos3 + 1])
                    if not np.any(np.isnan(dif)) and np.sum(dif > 0.3) > 5:
                        logRPCFed[startpos3:endpos3 + 1] = np.where(dif > 0.3, pcfed2, logRPCFed[startpos3:endpos3 + 1])
        logRPCFed = fillNA(logRPCFed, zeroIsNA=True)

        seg = rle_lengths(logRPCFed)

        # Initialize variables
        logRPCFed_new = np.array([], dtype=float)
        startprobe = 0
        prevlevel = 0

        for length in seg:
            endprobe = startprobe + length
            segment = logr[startprobe:endprobe]
            level = np.nanmean(segment)  # Compute mean ignoring NaNs

            # Making sure no NaN's get filled in
            if np.isnan(level):
                level = prevlevel
            else:
                prevlevel = level

            logRPCFed_new = np.concatenate((logRPCFed_new, np.full(length, level)))
            startprobe = endprobe

        logRPCFed = logRPCFed_new
        if len(np.unique(logRPCFed)) < 800:
            break

    bafPCFed = np.array(bafPCFed)
    Tumor_LogR_segmented = logRPCFed
    Tumor_BAF_segmented = 1 - bafPCFed

    tumor_logr_pcfed_output_dict = {key: Tumor_LogR_segmented[idx] for idx, key in enumerate(tumor_logr_dict.keys())}

    het = germline_genotypes_values == 'False'
    het_indices = np.where(het)[0]
    het_indices_keys = np.array(list(tumor_baf_dict.keys()))[het_indices]
    tumor_baf_pcfed_output_dict = {tuple(key): Tumor_BAF_segmented[idx] for idx, key in enumerate(het_indices_keys)}

    output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
    tumor_logr_pcfed_output = open(tumor_logr_pcfed_output_file, 'w')
    tumor_logr_pcfed_output.write(output_header)
    for key, value in tumor_logr_pcfed_output_dict.items():
        tumor_logr_pcfed_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
        tumor_logr_pcfed_output.write(tumor_logr_pcfed_string)
    tumor_baf_pcfed_output = open(tumor_baf_pcfed_output_file, 'w')
    tumor_baf_pcfed_output.write(output_header)
    for key, value in tumor_baf_pcfed_output_dict.items():
        tumor_baf_pcfed_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
        tumor_baf_pcfed_output.write(tumor_baf_pcfed_string)


def main():
    parser = ArgumentParser(description="Run ASPCF")

    parser.add_argument('--tumor_logr_file', type=str,
                        default=None,
                        help="Path of tumor sample LogR")

    parser.add_argument('--tumor_baf_file', type=str,
                        default=None,
                        help="Path of tumor sample BAF")

    parser.add_argument('--germline_genotypes_file', type=str,
                        default=None,
                        help="Path of germline genotypes")

    parser.add_argument('--tumor_logr_pcfed_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample PCFed LogR")

    parser.add_argument('--tumor_baf_pcfed_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample PCFed BAF")

    parser.add_argument('--penalty', type=int,
                        default=1000,
                        help="Penalty term")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")

    global args
    args = parser.parse_args()

    aspcf(args.tumor_logr_file, args.tumor_baf_file, args.germline_genotypes_file, args.tumor_logr_pcfed_output_file, args.tumor_baf_pcfed_output_file, args.penalty, args.sample_name)


if __name__ == "__main__":
    main()