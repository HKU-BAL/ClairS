# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import gc
import shlex
import tables
import numpy as np
import heapq
from random import random
from collections import defaultdict
from itertools import product

from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, IUPAC_base_to_num_dict as BASE2NUM

FILTERS = tables.Filters(complib='blosc:zstd', complevel=5)
shuffle_bin_size = 3000
PREFIX_CHAR_STR = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def setup_environment():
    gc.enable()


def batches_from(iterable, item_from, batch_size=1):
    iterable = iter(iterable)
    while True:
        chunk = []
        for _ in range(batch_size):
            try:
                chunk.append(item_from(next(iterable)))
            except StopIteration:
                yield chunk
                return
        yield chunk


def tensor_generator_from(tensor_file_path, batch_size, pileup, platform):
    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    tensor_shape = param.input_shape_dict[platform]
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        chrom, coord, seq, tensor, alt_info = row.split("\t")
        if pileup:
            tensor = np.array(tensor.split(), dtype=np.dtype(float_type))
            depth = int(alt_info.split('-', maxsplit=1)[0])
            max_depth = param.max_depth_dict[platform]
            # for extreme high coverage data, make sure we could have a truncated coverage
            if depth > 0 and depth > max_depth * 1.5:
                scale_factor = depth / max_depth
                tensor = tensor / scale_factor
        else:
            # need add padding if depth is lower than maximum depth.
            tensor = [int(item) for item in tensor.split()]
            tensor_depth = len(tensor) // tensor_shape[1] // tensor_shape[2]
            padding_depth = tensor_shape[0] - tensor_depth
            prefix_padding_depth = int(padding_depth / 2)
            suffix_padding_depth = padding_depth - int(padding_depth / 2)
            prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
            suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
            tensor = prefix_zero_padding + tensor + suffix_zero_padding
            tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = chrom + ":" + coord + ":" + seq
        return tensor, pos, seq, alt_info

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        alt_info_list = []
        for tensor, pos, seq, alt_info in batch:
            if seq[param.flankingBaseNum] not in BASE2NUM:
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            alt_info_list.append(alt_info)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= 0:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], alt_info_list[:current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


def variant_map_from(var_fn, tree, is_tree_empty):
    Y = {}
    miss_variant_set = set()
    if var_fn is None:
        return Y, miss_variant_set

    f = subprocess_popen(shlex.split("gzip -fdc %s" % (var_fn)))
    for row in f.stdout:
        columns = row.split()
        ctg_name, position_str = columns[0], columns[1]
        genotype1, genotype2 = columns[-2], columns[-1]
        key = ctg_name + ":" + position_str
        if genotype1 == '-1' or genotype2 == '-1':
            miss_variant_set.add(key)
            continue
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(position_str))):
            continue

        Y[key] = output_labels_from_vcf_columns(columns)

    f.stdout.close()
    f.wait()
    return Y, miss_variant_set


def write_table_dict(table_dict, normal_matrix, tumor_matrix, label, pos, total, normal_alt_info, tumor_alt_info,
                     tensor_shape, pileup, proportion=None, max_normal_depth=None, max_tumor_depth=None):
    normal_depth = len(normal_matrix) // tensor_shape[1] // tensor_shape[2]
    tumor_depth = len(tumor_matrix) // tensor_shape[1] // tensor_shape[2]
    tensor_depth = normal_depth + tumor_depth
    if max_normal_depth is not None and normal_depth > max_normal_depth:
        print("[WARNING] Skip adding high normal coverage data for position {}!".format('-'.join(pos.split(':')[:2])))
        return total
    if max_tumor_depth is not None and tumor_depth > max_tumor_depth:
        print("[WARNING] Skip adding high tumor coverage data for position {}!".format('-'.join(pos.split(':')[:2])))
        return total

    center_padding_depth = param.center_padding_depth
    padding_depth = tensor_shape[0] - tensor_depth - center_padding_depth
    prefix_padding_depth = int(padding_depth / 2)
    suffix_padding_depth = padding_depth - int(padding_depth / 2)
    prefix_zero_padding = ['0'] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
    center_zero_padding = ['0'] * center_padding_depth * tensor_shape[1] * tensor_shape[2]
    suffix_zero_padding = ['0'] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
    table_dict['input_matrix'].append(
        prefix_zero_padding + normal_matrix + center_zero_padding + tumor_matrix + suffix_zero_padding)
    table_dict['position'].append(pos)
    table_dict['label'].append(label)
    table_dict['normal_alt_info'].append(normal_alt_info)
    table_dict['tumor_alt_info'].append(tumor_alt_info)
    table_dict['proportion'].append(proportion)

    return total + 1


def update_table_dict():
    table_dict = {}
    table_dict['input_matrix'] = []
    table_dict['normal_alt_info'] = []
    table_dict['tumor_alt_info'] = []
    table_dict['position'] = []
    table_dict['label'] = []
    table_dict['proportion'] = []
    return table_dict


def write_table_file(table_file, table_dict, tensor_shape, label_size, float_type):
    """
    Write pileup or full alignment tensor into compressed bin file.
    table_dict: dictionary include all training information (tensor position, label, altnative bases).
    string: input tensor string, need add padding to meet the depth requirement.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set:  sometimes there will have true variant missing after downsampling reads.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    non_variant_subsample_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """

    try:
        input_matrix = np.array(table_dict['input_matrix'], np.dtype(float_type)).reshape([-1] + tensor_shape)
    except:
            return table_dict
    table_file.root.input_matrix.append(input_matrix)

    table_file.root.normal_alt_info.append(np.array(table_dict['normal_alt_info']).reshape(-1, 1))
    table_file.root.tumor_alt_info.append(np.array(table_dict['tumor_alt_info']).reshape(-1, 1))
    table_file.root.position.append(np.array(table_dict['position']).reshape(-1, 1))
    table_file.root.label.append(np.array(table_dict['label'], np.dtype(float_type)).reshape(-1, label_size))
    table_file.root.proportion.append(np.array(table_dict['proportion'], np.dtype('float32')).reshape(-1, 1))
    table_dict = update_table_dict()

    return table_dict


def print_bin_size(path, prefix=None):
    import tables
    import os
    total = 0
    for file_name in os.listdir(path):
        if prefix and not file_name.startswith(prefix):
            continue
        table = tables.open_file(os.path.join(path, file_name), 'r')
        print("[INFO] {} size is: {}".format(file_name, len(table.root.label)))
        total += len(table.root.label)
    print('[INFO] total: {}'.format(total))


def print_label(path):
    import tables
    import os
    total = 0
    total_germline = 0
    total_somatic = 0
    file_list = os.listdir(path) if not os.path.isfile(path) else [path]
    for file_name in file_list:
        table = tables.open_file(os.path.join(path, file_name), 'r')
        print("[INFO] {} size is: {}".format(file_name, len(table.root.label)))
        total += sum(table.root.label[:, 0])
        total_germline += sum(table.root.label[:, 1])
        total_somatic += sum(table.root.label[:, 2])
        table.close()
    print('[INFO] reference: {}, total germline:{}, total somatic:{}'.format(total, total_germline, total_somatic))


def check_bin_af_distrbution(path):
    import tables
    import os
    total = 0
    total_germline = 0
    total_somatic = 0
    normal_max = np.empty(shape=(0))
    tumor_max = np.empty(shape=(0))
    all_label = np.empty(shape=(0, 3))
    for fn_idx, file_name in enumerate(os.listdir(path)):
        table = tables.open_file(os.path.join(path, file_name), 'r')
        print("[INFO] {} size is: {}".format(file_name, len(table.root.label)))
        total += sum(table.root.label[:, 0])
        input_matrix_af_channel = table.root.input_matrix[:, :, 16, 4]
        normal_part = input_matrix_af_channel[:, :42]
        tumor_part = input_matrix_af_channel[:, 44:]
        label = table.root.label[:]
        normal_max = np.concatenate([normal_max, np.max(normal_part, axis=1)], axis=0)
        tumor_max = np.concatenate([tumor_max, np.max(tumor_part, axis=1)], axis=0)
        all_label = np.concatenate([all_label, label], axis=0)
        total_germline += sum(table.root.label[:, 1])
        total_germline += sum(table.root.label[:, 1])
        total_somatic += sum(table.root.label[:, 2])
        if fn_idx > 10:
            break
    if 0:
        import pandas as pd
        normal_max_df = pd.DataFrame(normal_max)
        normal = [True if np.argmax(x) == 0 else False for x in label]
        germline = [True if np.argmax(x) == 1 else False for x in label]
        somatic = [True if np.argmax(x) == 2 else False for x in label]
        input_matrix = np.array(table.root.input_matrix)
        position = np.array(table.root.position)
        normal_position = [p[0].decode().split(':')[0] for p in position if "ref" in p[0].decode()]
        germline_position = [p[0].decode().split(':')[0] for p in position if "homo_germline" in p[0].decode()]
        somatic_position = [p[0].decode().split(':')[0] for p in position if "homo_somatic" in p[0].decode()]
        normal_matrix = input_matrix[normal]
        germline_matrix = input_matrix[germline]
        somatic_matrix = input_matrix[somatic]

    print('[INFO] total: {}, total germline:{}, total somatic:{}'.format(total, total_germline, total_somatic))


def get_key_list(input_dict, shuffle=True):
    output_list = []
    for key, infos in input_dict.items():
        if 'normal' not in infos or 'tumor' not in infos:
            continue
        normal_index_list = range(len(infos['normal']))
        tumor_index_list = range(len(infos['tumor']))
        for x, y in product(normal_index_list, tumor_index_list):
            output_list.append((key, x, y))
    if shuffle == True:
        np.random.shuffle(output_list)
    return output_list


def bin_reader_generator_from(subprocess_process, Y, is_tree_empty, tree, miss_variant_set,
                              is_allow_duplicate_chr_pos=False, non_variant_subsample_ratio=1.0, is_tumor=False):
    """
    Bin reader generator for bin file generation.
    subprocess_list: a list includes all tensor generator of each tensor file.
    Y: dictionary (contig name: label information) to store all variant and non variant information.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set:  sometimes there will have true variant missing after downsampling reads.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    non_variant_subsample_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """
    dup_pos_end_flag = False
    pre_pos = None
    for row_idx, row in enumerate(subprocess_process.stdout):
        columns = row.split("\t")
        try:
            chrom, center_pos, seq, string, alt_info, tumor_tag, variant_type = columns
        except:
            continue
        # is_tumor = tumor_tag == 'tumor'
        alt_info = alt_info.rstrip()
        if not (is_tree_empty or is_region_in(tree, chrom, int(center_pos))):
            continue
        seq = seq.upper()
        if seq[param.flankingBaseNum] not in 'ACGT':
            continue
        key = chrom + ":" + center_pos
        is_reference = key not in Y

        if key in miss_variant_set:
            continue

        if is_reference and non_variant_subsample_ratio < 1.0 and random() >= non_variant_subsample_ratio:
            continue

        if pre_pos == center_pos:
            dup_pos_end_flag = False
        elif pre_pos and center_pos != pre_pos:
            dup_pos_end_flag = True
        pre_pos = center_pos
        key = chrom + ":" + center_pos
        yield (int(center_pos), key, is_tumor, string, alt_info, seq, variant_type, dup_pos_end_flag)

    subprocess_process.stdout.close()
    subprocess_process.wait()


def heapq_merge_generator_from(normal_bin_reader_generator, tumor_bin_reader_generator):
    tensor_infos_set = set()
    X = defaultdict(defaultdict)
    batch_count = 0
    normal_keys = set()
    tumor_keys = set()
    for tensor_infos in heapq.merge(normal_bin_reader_generator, tumor_bin_reader_generator):
        center_pos, key, is_tumor, string, alt_info, seq, somatic_flag, dup_pos_end_flag = tensor_infos
        if not is_tumor:
            normal_keys.add(key)
        else:
            tumor_keys.add(key)
        if key not in tensor_infos_set and is_tumor:
            continue
        tensor_infos_set.add(key)
        tumor_flag = 'tumor' if is_tumor else 'normal'
        tensor_list = string.split(" ")

        if batch_count >= shuffle_bin_size and is_tumor and dup_pos_end_flag:
            yield X, batch_count
            X = defaultdict(defaultdict)
            batch_count = 1

        if tumor_flag in X[key]:
            X[key][tumor_flag].append((tensor_list, alt_info, seq, somatic_flag))
        else:
            X[key][tumor_flag] = [(tensor_list, alt_info, seq, somatic_flag)]
        batch_count += 1

    if len(X):
        yield X, batch_count


def get_training_array(args,
                       normal_tensor_fn,
                       tumor_tensor_fn,
                       var_fn,
                       bed_fn,
                       bin_fn,
                       shuffle=True,
                       is_allow_duplicate_chr_pos=True,
                       chunk_id=None,
                       chunk_num=None,
                       platform='ont',
                       pileup=False,
                       maximum_non_variant_ratio=None,
                       phase_tumor=False,
                       candidate_details_fn_prefix=None,
                       merge_bins=False):
    tree = bed_tree_from(bed_file_path=bed_fn)
    is_tree_empty = len(tree.keys()) == 0
    Y, miss_variant_set = variant_map_from(var_fn, tree, is_tree_empty)

    global param
    float_type = 'int32'
    if pileup:
        import shared.param as param
    else:
        import shared.param as param
        float_type = 'int8'

    tensor_shape = param.input_shape_dict[platform]
    max_normal_depth = param.normal_matrix_depth_dict[platform]
    max_tumor_depth = param.tumor_matrix_depth_dict[platform]

    non_variant_subsample_ratio = maximum_non_variant_ratio if maximum_non_variant_ratio is not None else 1.0

    normal_tensor_list = []
    tumor_tensor_list = []
    if os.path.exists(tumor_tensor_fn) or os.path.exists(normal_tensor_fn):
        if not (os.path.exists(tumor_tensor_fn) and os.path.exists(normal_tensor_fn)):
            return 0
        normal_tensor_list.append(normal_tensor_fn)
        tumor_tensor_list.append(tumor_tensor_fn)
    else:
        tumor_tensor_info = tumor_tensor_fn.split('/')
        tumor_directry, file_prefix = '/'.join(tumor_tensor_info[:-1]), tumor_tensor_info[-1]

        normal_tensor_info = normal_tensor_fn.split('/')
        normal_directry, file_prefix = '/'.join(normal_tensor_info[:-1]), normal_tensor_info[-1]

        for file_name in os.listdir(tumor_directry):
            if file_name.startswith(file_prefix + '_') or file_name.startswith(
                    file_prefix + '.'):  # add '_.' to avoid add other prefix chr
                if os.path.exists(os.path.join(normal_directry, file_name)):
                    normal_tensor_list.append(os.path.join(normal_directry, file_name))
                    tumor_tensor_list.append(os.path.join(tumor_directry, file_name))

    tables.set_blosc_max_threads(64)
    int_atom = tables.Atom.from_dtype(np.dtype(float_type))
    float_atom = tables.Atom.from_dtype(np.dtype('float32'))
    string_atom = tables.StringAtom(itemsize=param.no_of_positions + 50)
    long_string_atom = tables.StringAtom(itemsize=30000)  # max alt_info length
    table_file = tables.open_file(bin_fn, mode='w', filters=FILTERS)
    table_file.create_earray(where='/', name='input_matrix', atom=float_atom, shape=[0] + tensor_shape,
                             filters=FILTERS)
    table_file.create_earray(where='/', name='position', atom=string_atom, shape=(0, 1), filters=FILTERS)
    table_file.create_earray(where='/', name='label', atom=float_atom, shape=(0, param.label_size), filters=FILTERS)
    table_file.create_earray(where='/', name='normal_alt_info', atom=long_string_atom, shape=(0, 1), filters=FILTERS)
    table_file.create_earray(where='/', name='tumor_alt_info', atom=long_string_atom, shape=(0, 1), filters=FILTERS)
    table_file.create_earray(where='/', name='proportion', atom=float_atom, shape=(0, 1), filters=FILTERS)

    table_dict = update_table_dict()
    total_compressed = 0
    total = 0
    for normal_tensor_fn, tumor_tensor_fn in zip(normal_tensor_list, tumor_tensor_list):

        normal_subprocess_process = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, normal_tensor_fn)))
        tumor_subprocess_process = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tumor_tensor_fn)))

        try:
            proportion = float(tumor_tensor_fn.split("_")[3])
        except:
            proportion = 1.0
        normal_bin_reader_generator = bin_reader_generator_from(subprocess_process=normal_subprocess_process,
                                                                Y=Y,
                                                                is_tree_empty=is_tree_empty,
                                                                tree=tree,
                                                                miss_variant_set=miss_variant_set,
                                                                is_allow_duplicate_chr_pos=is_allow_duplicate_chr_pos,
                                                                non_variant_subsample_ratio=non_variant_subsample_ratio)

        tumor_bin_reader_generator = bin_reader_generator_from(subprocess_process=tumor_subprocess_process,
                                                               Y=Y,
                                                               is_tree_empty=is_tree_empty,
                                                               tree=tree,
                                                               miss_variant_set=miss_variant_set,
                                                               is_allow_duplicate_chr_pos=is_allow_duplicate_chr_pos,
                                                               non_variant_subsample_ratio=non_variant_subsample_ratio,
                                                               is_tumor=True)

        total_compressed = 0
        total = 0
        for X, batch_total in heapq_merge_generator_from(normal_bin_reader_generator, tumor_bin_reader_generator):
            total += batch_total
            all_pos_list = get_key_list(X)

            for key, normal_index, tumor_index in all_pos_list:
                infos = X[key]
                normal_infos = infos['normal'][normal_index]
                tumor_infos = infos['tumor'][tumor_index]

                normal_tensor, nomral_alt_info, seq, variant_type = normal_infos
                variant_type = variant_type.rstrip()
                # use tumor alternative information instead of the normal one
                tumor_tensor, tumor_alt_info = tumor_infos[0], tumor_infos[1]
                pos_infos = key + ':' + seq + ':' + variant_type

                if variant_type == 'homo_somatic' or variant_type == "hetero_somatic":
                    label = [0, 0, 1]
                elif variant_type == 'homo_germline' or variant_type == 'hetero_germline':
                    label = [0, 1, 0]
                else:
                    label = [1, 0, 0]

                if args.use_reference_candidates_only and label[0] != 1:
                    continue

                total_compressed = write_table_dict(table_dict=table_dict,
                                                    normal_matrix=normal_tensor,
                                                    tumor_matrix=tumor_tensor,
                                                    label=label,
                                                    pos=pos_infos,
                                                    total=total_compressed,
                                                    normal_alt_info=nomral_alt_info,
                                                    tumor_alt_info=tumor_alt_info,
                                                    tensor_shape=tensor_shape,
                                                    pileup=pileup,
                                                    proportion=proportion,
                                                    max_normal_depth=max_normal_depth,
                                                    max_tumor_depth=max_tumor_depth
                                                    )

                if total_compressed % 500 == 0 and total_compressed > 0:
                    table_dict = write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)

        if total_compressed % 500 != 0 and total_compressed > 0:
            table_dict = write_table_file(table_file, table_dict, tensor_shape, param.label_size, float_type)

    table_file.close()
    print("[INFO] Compressed %d/%d tensor" % (total_compressed, total), file=sys.stderr)
