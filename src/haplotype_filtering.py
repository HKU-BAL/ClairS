import os
import subprocess
import concurrent.futures

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter
from shared.utils import str2bool, str_none

file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def get_base_list(columns):
    pileup_bases = columns[4]

    base_idx = 0
    base_list = []
    read_end_set = set()
    read_start_set = set()
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx].upper()
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance].upper()  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
            read_start_set.add(len(base_list) - 1)
        # skip $, the end of read
        if base == "$":
            read_end_set.add(len(base_list) - 1)
        base_idx += 1
    read_start_end_set = read_start_set if len(read_start_set) > len(read_end_set) else read_end_set
    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list, read_start_end_set


def extract_base(POS):
    pos = POS.pos
    ref_base = POS.reference_bases
    alt_base = POS.alternate_bases[0]
    tumor_bam_fn = args.tumor_bam_fn
    ref_fn = args.ref_fn
    samtools = args.samtools
    ctg_name = args.ctg_name if args.ctg_name is not None else POS.ctg_name
    min_mq = args.min_mq
    min_bq = args.min_bq

    pass_hetero = True
    pass_homo = True
    pass_hetero_both_side = True
    pass_read_start_end = True
    pass_bq = True
    flanking_list = POS.extra_infos
    hetero_germline_set = set()
    homo_germline_set = set()

    if not os.path.exists(tumor_bam_fn):
        tumor_bam_fn += ctg_name + '.bam'

    if flanking_list != "" and flanking_list is not None:
        hetero_germline_set, homo_germline_set = flanking_list

    flanking = args.flanking

    ctg_range = "{}:{}-{}".format(ctg_name, pos - flanking, pos + flanking + 1)
    samtools_command = "{} mpileup  --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-QNAME --output-extra HP ".format(
        samtools, min_mq, min_bq, ctg_range)

    tumor_samtools_command = samtools_command + tumor_bam_fn

    output = subprocess.run(tumor_samtools_command, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True).stdout.rstrip()

    from shared.utils import reference_sequence_from
    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools,
        fasta_file_path=ref_fn,
        regions=[ctg_range]
    )

    # tumor
    pos_dict = defaultdict(defaultdict)
    pos_counter_dict = defaultdict(defaultdict)
    hap_dict = defaultdict(int)
    all_read_start_end_set = set()

    homo_germline_pos_set = set([item[0] for item in homo_germline_set])
    hetero_germline_pos_set = set([item[0] for item in hetero_germline_set])

    for row in output.rstrip().split('\n'):
        columns = row.split('\t')
        read_name_list = columns[6].split(',')

        p = int(columns[1])
        ctg = columns[0]
        ctg = p if args.ctg_name is not None else (ctg, p)
        if ctg in hetero_germline_pos_set or p == pos:
            phasing_info = columns[7].split(',')
            for hap_idx, hap in enumerate(phasing_info):
                if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                    hap_dict[read_name_list[hap_idx]] = int(hap)

        base_counter, base_list, read_start_end_set = get_base_list(columns)

        # check union
        if len(read_start_end_set) >= len(base_list) * 0.2:
            all_read_start_end_set = all_read_start_end_set.union(
                set([read_name_list[r_idx] for r_idx in read_start_end_set]))

        pos_dict[p] = dict(zip(read_name_list, base_list))
        ref_base = reference_sequence[p - pos + flanking]

        if p == pos:
            bq_list = [ord(qual) - 33 for qual in columns[5]]
            alt_base_bq_set = [bq for key, value, bq in zip(read_name_list, base_list, bq_list) if
                               ''.join(value) == alt_base]
            # put into extract candidate
            if len(alt_base_bq_set) > 0 and sum(alt_base_bq_set) / len(alt_base_bq_set) <= 20 and float(POS.qual) < 0.9:
                pass_bq = False

            ALL_HAP_LIST = [0, 0, 0]
            HAP_LIST = [0, 0, 0]

            for rn in read_name_list:
                ALL_HAP_LIST[hap_dict[rn]] += 1

            alt_base_read_name_set = set(
                [key for key, value in zip(read_name_list, base_list) if ''.join(value) == alt_base])
            for rn in alt_base_read_name_set:
                HAP_LIST[hap_dict[rn]] += 1

        if len(base_counter) == 1 and base_counter[ref_base] > 0:
            continue
        pos_counter_dict[p] = base_counter

    # check del
    nor_del_base = 0
    del_base = len([key for key, value in pos_dict[pos].items() if ''.join(value) in '#*'])

    alt_base_read_name_set = set([key for key, value in pos_dict[pos].items() if ''.join(value) == alt_base])
    # near to read start end and have high overlap
    if len(all_read_start_end_set.intersection(alt_base_read_name_set)) >= 0.3 * len(alt_base_read_name_set):
        pass_read_start_end = False

    all_hap_counter = Counter([hap_dict[key] for key in pos_dict[pos].keys()])
    alt_hap_counter = Counter([hap_dict[key] for key in alt_base_read_name_set])

    all_hp0, all_hp1, all_hp2 = all_hap_counter[0], all_hap_counter[1], all_hap_counter[2]
    hp0, hp1, hp2 = alt_hap_counter[0], alt_hap_counter[1], alt_hap_counter[2]

    MAX = max(hp1, hp2)
    MIN = min(hp1, hp2)
    af = float(POS.af)
    if af < 0.1 and float(POS.qual) < 0.9:
        if hp1 * hp2 > 0 and MAX / MIN <= 10:
            pass_hetero_both_side = False

    is_phasable = hp1 * hp2 == 0 or (MAX / MIN >= 5 and (hp1 >= 3 or hp2 >= 3))
    hap_index = 0 if is_phasable else (1 if hp1 > hp2 else 2)
    # need to consider hap_0 hap_1 and hap_2

    # position with high overlap with current pos
    match_count = 0
    ins_length = 0
    all_base_dict = defaultdict(int)
    base_dict = defaultdict(int)
    alt_base_dict = defaultdict(int)

    for p, rb_dict in pos_dict.items():
        rb = reference_sequence[p - pos + flanking]
        read_alt_dict = pos_dict[p]
        if p == pos or p in homo_germline_pos_set or p in hetero_germline_pos_set:
            continue

        ins_length += sum(
            [min(len(v[1]) - 1, flanking * 2) for key, v in read_alt_dict.items() if len(v[1]) > 3 and v[1][0] == '+'])

        for k, v in read_alt_dict.items():
            hap = hap_dict[k]
            base = ''.join(v)
            all_base_dict[hap] += 1
            if base.upper() != rb:
                base_dict[hap] += 1
                if k in alt_base_read_name_set:
                    alt_base_dict[hap] += 1

        read_alt_dict = pos_dict[p]
        inter_set = set(read_alt_dict.keys()).intersection(alt_base_read_name_set)
        alt_list = []

        for key in inter_set:
            base = ''.join(read_alt_dict[key])
            if base != rb and base not in '#*':
                alt_list.append(base)
        alt_base_counter = sorted(Counter(alt_list).items(), key=lambda x: x[1], reverse=True)

        if len(alt_list) == 0 or (alt_base_counter[0][1] >= len(alt_base_read_name_set) * 1.2) or (
                alt_base_counter[0][1] <= len(alt_base_read_name_set) * 0.8):
            continue

        # cal all alt count in current position
        if pos_counter_dict[p][alt_base_counter[0][0]] >= alt_base_counter[0][1] * 1.3:
            continue
        match_count += 1

    for key, rb, ab in hetero_germline_set:
        read_alt_dict = pos_dict[p]
        # snp
        if len(rb) == 1 and len(ab) == 1:
            overlap_count = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
        # ins
        elif len(rb) == 1 and len(ab) > 1:
            overlap_count = set(
                [key for key, value in read_alt_dict.items() if ab[:2] in value[1][1:] and len(value[1]) > 1])
        # del
        elif len(rb) > 1 and len(ab) == 1:
            overlap_count = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])

        # in the same phased haplotype
        phased_overlap_count = set([key for key in overlap_count if hap_dict[key] == hap_index])
        if len(phased_overlap_count) == 0 or len(phased_overlap_count) / float(len(overlap_count)) < 0.5:
            continue

        inter_set = set([key for key in alt_base_read_name_set if hap_dict[key] == hap_index]).intersection(
            phased_overlap_count)
        if len(inter_set) == 0:
            pass_hetero = False
            break

    for p, rb, ab in homo_germline_set:
        read_alt_dict = pos_dict[p]

        # is the homo confident
        if len(rb) == 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
        elif len(rb) == 1 and len(ab) > 1:
            homo_alt_key = set(
                [key for key, value in read_alt_dict.items() if (ab[1:2] in value[1][1:]) and len(value[1]) > 1])
        elif len(rb) > 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])

        homo_hap_counter = Counter([hap_dict[key] for key in homo_alt_key])
        all_homo_hap_counter = Counter([hap_dict[key] for key in read_alt_dict.keys()])

        all_hp0, all_hp1, all_hp2 = all_homo_hap_counter[0], all_homo_hap_counter[1], all_homo_hap_counter[2]
        hp0, hp1, hp2 = homo_hap_counter[0], homo_hap_counter[1], homo_hap_counter[2]
        af = (hp0 + hp1 + hp2) / float(all_hp0 + all_hp1 + all_hp2)

        def phasble(all_hap_list, hap_list):
            all_hp0, all_hp1, all_hp2 = all_hap_list
            if all_hp1 * all_hp2 == 0:
                return False
            hp0, hp1, hp2 = hap_list
            MAX = max(hp1, hp2)
            MIN = min(hp1, hp2)
            if hp1 * hp2 > 0 and MAX / MIN <= 10:
                return False
            return True

        is_homo_alt_phasable = phasble([all_hp0, all_hp1, all_hp2], [hp0, hp1, hp2])

        if af < 0.75 or is_homo_alt_phasable:
            continue

        inter_set = set(list(read_alt_dict.keys())).intersection(alt_base_read_name_set)
        if len(inter_set) == 0:
            continue
        if len(rb) == 1 and len(ab) == 1:
            overlap_count = set(
                [key for key, value in read_alt_dict.items() if ''.join(value) == ab and key in inter_set])
        elif len(rb) == 1 and len(ab) > 1:
            overlap_count = set([key for key, value in read_alt_dict.items() if
                                 (ab[1:2] in value[1][1:]) and len(value[1]) > 1 and key in inter_set])
        elif len(rb) > 1 and len(ab) == 1:
            overlap_count = set(
                [key for key, value in read_alt_dict.items() if '-' in ''.join(value) and key in inter_set])
        if len(overlap_count) == 0 or len(overlap_count) / len(inter_set) < 0.2:
            pass_homo = False
            break

    key = (ctg_name, pos) if args.ctg_name is None else pos
    return key, pass_hetero, pass_homo, pass_hetero_both_side, match_count, ins_length, pass_read_start_end, alt_base_read_name_set, (
        nor_del_base, del_base, all_base_dict, base_dict, alt_base_dict, pass_bq, ALL_HAP_LIST, HAP_LIST)


def realign_variants(args):
    ctg_name = args.ctg_name
    threads = args.threads
    apply_post_processing = args.apply_post_processing
    pileup_vcf_fn = args.pileup_vcf_fn
    fa_input_vcf_fn = args.full_alignment_vcf_fn
    germline_vcf_fn = args.germline_vcf_fn
    flanking = args.flanking
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        subprocess.run("mkdir -p {}".format(output_dir), shell=True)

    pileup_output_vcf_fn = os.path.join(output_dir, "pileup_filter.vcf")
    fa_output_vcf_fn = os.path.join(output_dir, "full_alignment_filter.vcf")
    if not apply_post_processing:
        subprocess.run("ln -sf {} {}".format(pileup_vcf_fn, pileup_output_vcf_fn), shell=True)
        subprocess.run("ln -sf {} {}".format(fa_input_vcf_fn, fa_output_vcf_fn), shell=True)
        return

    germine_input_vcf_reader = VcfReader(vcf_fn=germline_vcf_fn, ctg_name=ctg_name, show_ref=False, keep_row_str=False,
                                         filter_tag="PASS", save_header=False,
                                         skip_genotype=False)
    germine_input_vcf_reader.read_vcf()
    germline_input_variant_dict = germine_input_vcf_reader.variant_dict

    germline_gt_list = []
    # only keep homo germline
    for key in list(germline_input_variant_dict.keys()):
        if sum(germline_input_variant_dict[key].genotype) == 1:
            germline_gt_list.append((key, 1))
        elif sum(germline_input_variant_dict[key].genotype) == 2:
            germline_gt_list.append((key, 2))

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn, ctg_name=ctg_name, show_ref=False, keep_row_str=True,
                                 discard_indel=True,
                                 filter_tag=args.input_filter_tag, save_header=True,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    input_vcf_reader = VcfReader(vcf_fn=fa_input_vcf_fn, ctg_name=ctg_name, show_ref=False, keep_row_str=True,
                                 discard_indel=True,
                                 filter_tag=args.input_filter_tag, save_header=True,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    fa_variant_dict = input_vcf_reader.variant_dict

    input_variant_dict = defaultdict()
    for k, v in fa_variant_dict.items():
        if v.filter != "PASS":
            continue
        if len(pileup_variant_dict) and (k not in pileup_variant_dict or pileup_variant_dict[k].filter != "PASS"):
            continue
        input_variant_dict[k] = v

    print("Total Input: ", ctg_name, len(pileup_variant_dict), len(fa_variant_dict), len(input_variant_dict))

    p_vcf_writer = VcfWriter(vcf_fn=pileup_output_vcf_fn, ctg_name=ctg_name, show_ref_calls=True)
    f_vcf_writer = VcfWriter(vcf_fn=fa_output_vcf_fn, ctg_name=ctg_name, show_ref_calls=True)

    for key, POS in input_variant_dict.items():
        pos = key if args.ctg_name is not None else key[1]
        hetero_flanking_list = []
        homo_flanking_list = []
        for gk, gt in germline_gt_list:
            p = gk if args.ctg_name is not None else gk[1]
            if p > pos + flanking:
                break
            if p > pos - flanking and p != pos:
                ref_base = germline_input_variant_dict[gk].reference_bases
                alt_base = germline_input_variant_dict[gk].alternate_bases[0]
                if gt == 1:
                    hetero_flanking_list.append((p, ref_base, alt_base))
                else:
                    homo_flanking_list.append((p, ref_base, alt_base))
        POS.extra_infos = [set(hetero_flanking_list), set(homo_flanking_list)]

    co_exist_fail_pos_set = set()
    complex_indel_fail_pos_set = set()
    fail_hetero_set = set()
    fail_homo_set = set()
    fail_hetero_both_side_set = set()
    fail_pass_read_start_end_set = set()
    fail_bq_set = set()
    phase_dict = defaultdict()
    fail_rn_set = set()

    ctg_name_list = [ctg_name] if args.ctg_name is not None else set([item[0] for item in input_variant_dict.keys()])
    for ctg in ctg_name_list:
        POS_list = list([v for k, v in input_variant_dict.items() if
                         k[0] == ctg]) if args.ctg_name is None else input_variant_dict.values()

        rn_dict = defaultdict(list)

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as exec:
            for result in exec.map(extract_base, POS_list):
                key, pass_hetero, pass_homo, pass_hetero_both_side, match_count, indel_length, pass_read_start_end, alt_base_read_name_set, extra = result
                ctg_name = key[0] if args.ctg_name is None else args.ctg_name
                pos = key[1] if args.ctg_name is None else key
                # rn_dict[(ctg_name, int(pos))] = alt_base_read_name_set
                if match_count > 2:
                    co_exist_fail_pos_set.add((ctg_name, pos))

                nor_del_base, del_base, all_base_dict, base_dict, alt_base_dict, pass_bq, ALL_HAP_LIST, HAP_LIST = extra

                if indel_length > 300 and float(input_variant_dict[key].qual) < 0.9:
                    complex_indel_fail_pos_set.add((ctg_name, pos))

                if pass_hetero is False:
                    fail_hetero_set.add((ctg_name, pos))
                if pass_homo is False:
                    fail_homo_set.add((ctg_name, pos))

                if pass_hetero_both_side is False:
                    fail_hetero_both_side_set.add((ctg_name, pos))

                if pass_read_start_end is False:
                    fail_pass_read_start_end_set.add((ctg_name, pos))

                if pass_bq is False:
                    fail_bq_set.add((ctg_name, pos))

                phase_dict[(ctg_name, pos)] = ALL_HAP_LIST + HAP_LIST

        max_distance = 100000
        key_list = sorted(rn_dict.keys(), key=lambda x: (x[0], int(x[1])))

        for idx, k in enumerate(key_list):
            ctg, pos = k
            rn = rn_dict[k]
            count = 0
            for i in range(len(key_list)):
                if i == idx:
                    continue
                tmp_key = key_list[i]
                c, p = tmp_key
                if c != ctg:
                    break
                if int(pos) - int(p) > max_distance:
                    continue

                if int(p) - int(pos) > max_distance:
                    break

                tmp_rn = rn_dict[tmp_key]
                if len(rn.intersection(tmp_rn)) >= len(rn) * 0.5:
                    count += 1
            if count > 0:
                fail_rn_set.add(k)

    for key in sorted(fa_variant_dict.keys()):
        row_str = fa_variant_dict[key].row_str.rstrip()
        ctg_name = key[0] if args.ctg_name is None else args.ctg_name
        pos = key[1] if args.ctg_name is None else key

        if (ctg_name, pos) in co_exist_fail_pos_set:
            row_str = row_str.replace("PASS", "Low_co_exist")
            row_str = row_str + ':' + "Low_co_exist"

        if (ctg_name, pos) in complex_indel_fail_pos_set:
            row_str = row_str.replace("PASS", "Low_complex_indel")
            row_str = row_str + ':' + "Low_complex_indel"

        if (ctg_name, pos) in fail_hetero_set:
            row_str = row_str.replace("PASS", "Low_hetero")
            row_str = row_str + ':' + "Low_hetero"

        if (ctg_name, pos) in fail_homo_set:
            row_str = row_str.replace("PASS", "Low_homo")
            row_str = row_str + ':' + "Low_homo"

        if (ctg_name, pos) in fail_hetero_both_side_set:
            row_str = row_str.replace("PASS", "Low_hetero_both")
            row_str = row_str + ':' + "Low_hetero_both"

        if (ctg_name, pos) in fail_pass_read_start_end_set:
            row_str = row_str + ':' + "Low_read_start_end"

        if (ctg_name, pos) in fail_bq_set:
            row_str = row_str.replace("PASS", "Low_bq")
            row_str = row_str + ':' + "Low_bq"

        if (ctg_name, pos) in fail_bq_set:
            row_str = row_str.replace("PASS", "Low_distance")
            row_str = row_str + ':' + "Low_distance"

        phasing_info = ' '.join([str(item) for item in phase_dict[(ctg_name, pos)]]) if (ctg_name,
                                                                                         pos) in phase_dict else "0 0 0 0 0 0"
        row_str += "\t" + phasing_info
        f_vcf_writer.vcf_writer.write(row_str + '\n')

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        ctg_name = key[0] if args.ctg_name is None else args.ctg_name
        pos = key[1] if args.ctg_name is None else key

        if (ctg_name, pos) in co_exist_fail_pos_set:
            row_str = row_str.replace("PASS", "Low_co_exist")
            row_str = row_str + ':' + "Low_co_exist"

        if (ctg_name, pos) in complex_indel_fail_pos_set:
            row_str = row_str.replace("PASS", "Low_complex_indel")
            row_str = row_str + ':' + "Low_complex_indel"

        if (ctg_name, pos) in fail_hetero_set:
            row_str = row_str.replace("PASS", "Low_hetero")
            row_str = row_str + ':' + "Low_hetero"

        if (ctg_name, pos) in fail_homo_set:
            row_str = row_str.replace("PASS", "Low_homo")
            row_str = row_str + ':' + "Low_homo"

        if (ctg_name, pos) in fail_hetero_both_side_set:
            row_str = row_str.replace("PASS", "Low_hetero_both")
            row_str = row_str + ':' + "Low_hetero_both"

        if (ctg_name, pos) in fail_pass_read_start_end_set:
            row_str = row_str + ':' + "Low_read_start_end"

        if (ctg_name, pos) in fail_bq_set:
            row_str = row_str.replace("PASS", "Low_bq")
            row_str = row_str + ':' + "Low_bq"

        if (ctg_name, pos) in fail_bq_set:
            row_str = row_str.replace("PASS", "Low_distance")
            row_str = row_str + ':' + "Low_distance"

        phasing_info = ' '.join([str(item) for item in phase_dict[(ctg_name, pos)]]) if (ctg_name,
                                                                                         pos) in phase_dict else "0 0 0 0 0 0"
        row_str += "\t" + phasing_info
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()
    f_vcf_writer.close()


def main():
    parser = ArgumentParser(description="Haplotype filter for long read data")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--normal_bam_fn', type=str, default=None,
                        help="Sorted normal BAM file input")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Full-alignment VCF input")

    parser.add_argument('--germline_vcf_fn', type=str, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output vcf directory")

    parser.add_argument('--apply_post_processing', type=str2bool, default=True,
                        help="Apply post processing in calling")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")

    parser.add_argument('--threads', type=int, default=4,
                        help="Max #threads to be used")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help='VCF FILTER tag for input VCF')

    # options for advanced users
    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    ## test using one position
    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    ## flakning window size to process
    parser.add_argument('--flanking', type=int, default=100,
                        help=SUPPRESS)

    global args
    args = parser.parse_args()

    realign_variants(args)


if __name__ == "__main__":
    main()