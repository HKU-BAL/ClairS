import logging
import tables
import sys

import os
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from time import time
from subprocess import run
from argparse import ArgumentParser, SUPPRESS
from torch.optim.lr_scheduler import StepLR
from torch.utils.tensorboard import SummaryWriter
from torchsummary import summary
from tqdm import tqdm
from collections import defaultdict

from shared.utils import str2bool
import shared.param as param
import clair_somatic.model as model_path
from shared.vcf import VcfReader

# reuqired package  torchsummary, tqdm tables,  einops
logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '512'
os.environ['NUMEXPR_NUM_THREADS'] = '256'
batch_size = 64
gamma = 0.7
seed = 100
random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True


class FocalLoss(nn.Module):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, alpha=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma

    def forward(self, input, target):
        y_pred, y_true = input, target
        y_pred = torch.clamp(y_pred, min=1e-9, max=1 - 1e-9)
        cross_entropy = -y_true * torch.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight

        reduce_fl = torch.mean(FCLoss)
        return reduce_fl

def cal_metrics(tp, fp, fn):

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return round(precision, 6), round(recall, 6), round(f1_score, 6)


def get_label_task(label, label_shape_cum, task):
    if task == 0:
        return label[:label_shape_cum[task]]
    elif task == len(label_shape_cum) - 1:
        return label[label_shape_cum[task - 1]:]
    else:
        return label[label_shape_cum[task - 1]:label_shape_cum[task]]


def cal_class_weight(samples_per_cls, no_of_classes, beta=0.999):
    effective_num = 1.0 - np.power(beta, samples_per_cls)
    cls_weights = (1.0 - beta) / np.array(effective_num)
    cls_weights = cls_weights / np.sum(cls_weights) * no_of_classes
    return cls_weights


def get_chunk_list(chunk_offset, train_chunk_num, chunks_per_batch=10, training_dataset_percentage=None):
    """
    get chunk list for training and validation data. we will randomly split training and validation dataset,
    all training data is directly acquired from various tensor bin files.

    """
    need_split_validation_data = training_dataset_percentage is not None
    all_shuffle_chunk_list = []
    training_chunk_list, validation_chunk_list = [], []
    for bin_idx, chunk_num in enumerate(chunk_offset):
        current_chunk_list = [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
        all_shuffle_chunk_list += current_chunk_list
        if need_split_validation_data:
            buffer_chunk_num = chunks_per_batch
            if chunk_num < buffer_chunk_num:
                training_chunk_list += [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
                continue

            training_chunk_num = int((chunk_num - buffer_chunk_num) * training_dataset_percentage)
            validation_chunk_num = int(chunk_num - buffer_chunk_num - training_chunk_num)
            if training_chunk_num > 0:
                training_chunk_list += current_chunk_list[:training_chunk_num]
            if validation_chunk_num > 0:
                validation_chunk_list += current_chunk_list[-validation_chunk_num:]

    if need_split_validation_data:
        return np.array(training_chunk_list), np.array(validation_chunk_list)

    return np.array(all_shuffle_chunk_list[:train_chunk_num]), np.array(all_shuffle_chunk_list[train_chunk_num:])


def exist_file_prefix(exclude_training_samples, f):
    for prefix in exclude_training_samples:
        if prefix in f:
            return True
    return False


def pass_chr(fn, ctg_name_list):
    if ctg_name_list is None or len(ctg_name_list) == 0:
        return True
    in_testing_chr = False
    for ctg_name in ctg_name_list:
        if '_' + ctg_name + '.' in fn:
            return True
    return False

def cal_metrics(tp, fp, fn):

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return round(precision, 6), round(recall, 6), round(f1_score, 6)

def find_max_candidates(pos, alt_info):
    depth, alt_str = alt_info.split('-')[:2]
    # depth, alt_str, af_str = alt_info.split('-')
    alt_list = alt_str.split(' ')
    alt_list = list(zip(alt_list[::2], [int(item) for item in alt_list[1::2]])) if len(alt_list) else []
    alt_list = sorted(alt_list, key = lambda x: x[1], reverse=True)
    if len(alt_list) == 0:
        return False, False, False, None
    best_match_alt = alt_list[0][0]
    is_snp = best_match_alt[0] == "X"
    is_ins = best_match_alt[0] == "I"
    is_del = best_match_alt[0] == "D"
    return is_snp, is_ins, is_del, best_match_alt

def get_max_af(input_matrix):

    normal_tumor_channel = input_matrix[:,6,:]
    is_normal = normal_tumor_channel == 0.5
    is_tumor = normal_tumor_channel == 1.0
    max_normal = np.max(input_matrix[:,4,:] * is_normal, axis=1)
    max_tumor = np.max(input_matrix[:,4,:] * is_tumor, axis=1)
    # all_max_normal = []
    # all_max_tumor = []
    # for item in input_matrix:
    #     normal_tumor_channel = item[6,:]
    #     is_normal = [True if x == 0.5 else False for x in normal_tumor_channel]
    #     is_tumor = [True if x == 1.0 else False for x in normal_tumor_channel]
    #     max_normal = np.max(item[4, :][is_normal])
    #     max_tumor = np.max(item[4,:][is_tumor])
    #     all_max_normal.append(max_normal)
    #     all_max_tumor.append(max_tumor)
    return max_normal, max_tumor

def predict_model(args):
    apply_focal_loss = False
    use_resnet = args.use_resnet
    platform = args.platform
    output_dir = args.output_dir
    ctg_name_string = args.ctgName
    chkpnt_fn = args.chkpnt_fn
    ochk_prefix = args.ochk_prefix
    add_writer = args.add_writer
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num

    add_validation_dataset = args.random_validation or (args.validation_fn is not None)
    validation_fn = args.validation_fn
    ctg_name_list = ctg_name_string.split(',') if ctg_name_string is not None else []
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    label_size, label_shape = param.label_size, param.label_shape
    test_all_pos = args.test_all_pos
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""
    if ochk_prefix and not os.path.exists(ochk_prefix):
        output = run('mkdir {}'.format(ochk_prefix), shell=True)
        print("[INFO] Model path empty, create folder")

    if output_dir and not os.path.exists(output_dir):
        output = run('mkdir {}'.format(output_dir), shell=True)
        print("[INFO] Output path empty, create folder")

    unified_vcf_fn = args.unified_vcf_fn
    unified_vcf_reader_dict = defaultdict()
    for ctg_name in ctg_name_list:
        unified_vcf_reader = VcfReader(vcf_fn=unified_vcf_fn + ctg_name, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        unified_vcf_reader_dict[ctg_name] = unified_vcf_reader.variant_dict

    if add_writer:
        writer = SummaryWriter("{}/log".format(ochk_prefix))

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = model_path.CvT(
        num_classes=3,
        s1_emb_dim=16,  # stage 1 - dimension
        s1_emb_kernel=3,  # stage 1 - conv kernel
        s1_emb_stride=2,  # stage 1 - conv stride
        s1_proj_kernel=3,  # stage 1 - attention ds-conv kernel size
        s1_kv_proj_stride=2,  # stage 1 - attention key / value projection stride
        s1_heads=1,  # stage 1 - heads
        s1_depth=1,  # stage 1 - depth
        s1_mlp_mult=4,  # stage 1 - feedforward expansion factor
        s2_emb_dim=64,  # stage 2 - (same as above)
        s2_emb_kernel=3,
        s2_emb_stride=2,
        s2_proj_kernel=3,
        s2_kv_proj_stride=2,
        s2_heads=3,
        s2_depth=2,
        s2_mlp_mult=4,
        s3_emb_dim=128,  # stage 3 - (same as above)
        s3_emb_kernel=3,
        s3_emb_stride=2,
        s3_proj_kernel=3,
        s3_kv_proj_stride=2,
        s3_heads=4,
        s3_depth=3,
        s3_mlp_mult=4,
        dropout=0.,
        depth=param.max_depth,
        width=param.no_of_positions,
        dim=param.channel_size,
        apply_softmax=True if apply_focal_loss else False
    ).to(device)

    if use_resnet:
        model = model_path.ResNet().to(device)
    if chkpnt_fn is not None:
        model = torch.load(chkpnt_fn)
    input = torch.ones(size=(100, param.channel_size, param.max_depth, param.no_of_positions)).to(device)
    output = model(input)
    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape

    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    assert batch_size % chunk_size == 0
    chunks_per_batch = batch_size // chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    bin_list = os.listdir(args.bin_fn)
    bin_list = [f for f in bin_list if pass_chr(f, ctg_name_list) and not exist_file_prefix(exclude_training_samples, f)]
    bin_list = sorted(bin_list)
    total_bins = len(bin_list)
    if chunk_id is not None and chunk_num is not None:
        chunk_bin_size = total_bins // chunk_num if total_bins % chunk_num == 0 else total_bins // chunk_num + 1
        chunk_start_pos = chunk_id * chunk_bin_size
        chunk_end_pos = chunk_start_pos + chunk_bin_size
        bin_list = bin_list[chunk_start_pos:chunk_end_pos]
    # bin_list = bin_list[:3]
    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))

    def populate_dataset_table(file_list, file_path):
        chunk_offset = np.zeros(len(file_list), dtype=int)
        table_dataset_list = []
        for bin_idx, bin_file in enumerate(file_list):
            table_dataset = tables.open_file(os.path.join(file_path, bin_file), 'r')
            table_dataset_list.append(table_dataset)
            chunk_num = (len(table_dataset.root.label) - batch_size) // chunk_size
            chunk_offset[bin_idx] = chunk_num
        return table_dataset_list, chunk_offset

    table_dataset_list, chunk_offset = populate_dataset_table(bin_list, args.bin_fn)

    validate_table_dataset_list = []
    if validation_fn:
        val_list = os.listdir(validation_fn)
        logging.info("[INFO] total {} validation bin files: {}".format(len(val_list), ','.join(val_list)))
        validate_table_dataset_list, validate_chunk_offset = populate_dataset_table(val_list, args.validation_fn)

        train_chunk_num = int(sum(chunk_offset))
        train_shuffle_chunk_list, _ = get_chunk_list(chunk_offset, train_chunk_num)

        validate_chunk_num = int(sum(validate_chunk_offset))
        validate_shuffle_chunk_list, _ = get_chunk_list(validate_chunk_offset, validate_chunk_num)
        total_chunks = train_chunk_num + validate_chunk_num
    else:
        total_chunks = int(sum(chunk_offset))
        training_dataset_percentage = param.trainingDatasetPercentage if add_validation_dataset else None
        if add_validation_dataset:
            total_batches = total_chunks // chunks_per_batch
            validate_chunk_num = int(max(1., np.floor(total_batches * (1 - training_dataset_percentage))) * chunks_per_batch)
            # +++++++++++++**----
            # +:training *:buffer -:validation
            # distribute one batch data as buffer for each bin file, avoiding shifting training data to validation data
            train_chunk_num = int(total_chunks - validate_chunk_num)
        else:
            train_chunk_num = total_chunks
        train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_chunk_num, chunks_per_batch, training_dataset_percentage)
        train_chunk_num = len(train_shuffle_chunk_list)
        validate_chunk_num = len(validate_shuffle_chunk_list)
    train_data_size = train_chunk_num * chunk_size
    validate_data_size = validate_chunk_num * chunk_size

    batch_size = param.predictBatchSize
    def DataGenerator_all(x):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """

        for dataset in x:
            data_size = len(dataset.root.label)
            batch_num = data_size // batch_size if data_size % batch_size == 0 else data_size // batch_size + 1

            for batch_idx in range(batch_num):
                current_batch_size = batch_size if batch_idx != batch_num - 1 else data_size - batch_size * batch_idx
                input_matrix = np.array(dataset.root.input_matrix[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                label = np.array(dataset.root.label[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                positions = np.array(dataset.root.position[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                alt_infos = np.array(dataset.root.alt_info[batch_idx * batch_size:(batch_idx + 1) * batch_size])

                contig_infos = dataset.filename.split('/')[-1].split('_')[1].split('.')[0]
                ctg_name_infos = np.array([[contig_infos]] * current_batch_size)

                input_tensor = torch.from_numpy(np.transpose(input_matrix.astype(float), (0,3,1,2))/100.0).to(device)
                # tumor_tensor = torch.from_numpy(tumor_matrix)
                label_tensor = label
                yield input_tensor, label_tensor, positions, alt_infos, ctg_name_infos
    validate_dataset = DataGenerator_all(table_dataset_list)

    print("[INFO] The size of dataset: {}".format(train_data_size))
    print("[INFO] The training batch size: {}".format(batch_size))
    print("[INFO] The training learning_rate: {}".format(learning_rate))
    print("[INFO] The output model folder: {}".format(ochk_prefix))
    if chkpnt_fn is not None:
        model = torch.load(chkpnt_fn)

    val_total = 0
    predictions = np.empty(shape=(0, label_size))
    normal_logits = np.empty(shape=(0, 64))
    tumor_logits = np.empty(shape=(0, 64))
    labels = np.empty(shape=(0, label_size))
    predict_germline = False
    prefix = 'germline_' if predict_germline else ""
    if chunk_id is not None: prefix  += str(chunk_id)
    fp_vcf = open(os.path.join(output_dir, prefix +'fp.vcf'), 'w')
    fn_vcf = open(os.path.join(output_dir, prefix+'fn.vcf'), 'w')
    tp_vcf = open(os.path.join(output_dir, prefix+'tp.vcf'), 'w')

    tp_snp, tp_ins, tp_del, fp_snp, fp_ins, fp_del, fn_snp, fn_ins, fn_del, fp_snp_truth, fp_ins_truth, fp_del_truth = 0,0,0,0,0,0,0,0,0,0,0,0
    tp_count = 0
    for val_data, val_label, val_positions, val_alt_infos, val_ctg_name_infos in validate_dataset:
        with torch.no_grad():
            prediction = model(val_data.float())
            prediction = prediction.cpu().numpy()
        val_total += len(prediction)
        # predictions = np.concatenate([predictions, prediction], axis = 0)
        y_pred = np.argmax(prediction, axis=1)
        y_truth = np.argmax(val_label, axis=1)
        #np where to get index
        arg_index = 1 if predict_germline else 2
        fp = [True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)]
        fn = [True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)]
        tp = [True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)]

        input_matrix = np.max(val_data.cpu().numpy(), axis=3)
        max_normal_af, max_tumor_af = get_max_af(input_matrix)
        try:
            fp_pos_array = [item[0].numpy().decode().split(':') for item in val_positions[fp]]
            fn_pos_array = [item[0].numpy().decode().split(':') for item in val_positions[fn]]
            tp_pos_array = [item[0].numpy().decode().split(':') for item in val_positions[tp]]

            fp_chr_array = [item[0].numpy().decode() for item in val_ctg_name_infos[fp]]
            fn_chr_array = [item[0].numpy().decode() for item in val_ctg_name_infos[fn]]
            tp_chr_array = [item[0].numpy().decode() for item in val_ctg_name_infos[tp]]

            fp_alt_array = [item[0].numpy().decode() for item in val_alt_infos[fp]]
            fn_alt_array = [item[0].numpy().decode() for item in val_alt_infos[fn]]
            tp_alt_array = [item[0].numpy().decode() for item in val_alt_infos[tp]]
        except:
            fp_pos_array = [item[0].decode().split(':') for item in val_positions[fp]]
            fn_pos_array = [item[0].decode().split(':') for item in val_positions[fn]]
            tp_pos_array = [item[0].decode().split(':') for item in val_positions[tp]]

            fp_chr_array = [item[0] for item in val_ctg_name_infos[fp]]
            fn_chr_array = [item[0] for item in val_ctg_name_infos[fn]]
            tp_chr_array = [item[0] for item in val_ctg_name_infos[tp]]

            fp_alt_array = [item[0].decode() for item in val_alt_infos[fp]]
            fn_alt_array = [item[0].decode() for item in val_alt_infos[fn]]
            tp_alt_array = [item[0].decode() for item in val_alt_infos[tp]]

        fp_normal_af_array = max_normal_af[fp]
        fn_normal_af_array = max_normal_af[fn]
        tp_normal_af_array = max_normal_af[tp]
        fp_tumor_af_array = max_tumor_af[fp]
        fn_tumor_af_array = max_tumor_af[fn]
        tp_tumor_af_array = max_tumor_af[tp]

        fp_pos_array = [(item[0],item[-1]) for item in fp_pos_array]
        fn_pos_array = [(item[0],item[-1]) for item in fn_pos_array]
        tp_pos_array = [(item[0],item[-1]) for item in tp_pos_array]

        all_pos_array = [fp_pos_array, fn_pos_array, tp_pos_array]
        all_chr_array = [fp_chr_array, fn_chr_array, tp_chr_array]
        all_alt_array = [fp_alt_array, fn_alt_array, tp_alt_array]
        all_normal_af_array = [fp_normal_af_array, fn_normal_af_array, tp_normal_af_array]
        all_tumor_af_array = [fp_tumor_af_array, fn_tumor_af_array, tp_tumor_af_array]

        all_vcf_fn = [fp_vcf, fn_vcf, tp_vcf]
        candidate_types = ['fp', 'fn', 'tp']

        for pos_array, chr_array, vcf_fn, alt_array, ct, normal_af_array, tumor_af_array in zip(all_pos_array, all_chr_array, all_vcf_fn, all_alt_array, candidate_types, all_normal_af_array, all_tumor_af_array):
            # print(len(pos_array), len(chr_array), len(alt_array), len(ct))
            for (pos, variant_type), contig, alt_info, n_af, t_af in zip(pos_array, chr_array, alt_array, normal_af_array, tumor_af_array):
                vcf_format = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:NAF:TAF:VT\t%s:%d:%d:%.4f:%.4f:%s" % (
                    contig,
                    int(pos),
                    "A",
                    "A",
                    10,
                    'PASS',
                    '.',
                    "0/0",
                    10,
                    10,
                    n_af,
                    t_af,
                    variant_type)
                vcf_fn.write(vcf_format + '\n')

                # calculate indel and snp information
                is_snp, is_ins, is_del, best_match_alt = find_max_candidates(pos, alt_info)

                if contig in unified_vcf_reader_dict and int(pos) in unified_vcf_reader_dict[contig]:
                    ref_base = unified_vcf_reader_dict[contig][int(pos)].reference_bases
                    alt_base = unified_vcf_reader_dict[contig][int(pos)].alternate_bases[0]
                    is_snp_truth = len(ref_base) == 1 and len(alt_base) == 1
                    is_ins_truth = len(ref_base) < len(alt_base)
                    is_del_truth = len(ref_base) > len(alt_base)
                    tp_count += 1

                # if is_snp_truth and not is_snp:
                #     print (pos, ref_base, alt_base,best_match_alt, alt_info)

                # if best_match_alt is None:
                #     is_snp = is_snp_truth
                #     is_ins = is_ins_truth
                #     is_del = is_del_truth

                    # print(ref_base, alt_base,best_match_alt)

                tp_snp = tp_snp + 1 if ct == 'tp' and is_snp and is_snp_truth else tp_snp
                tp_ins = tp_ins + 1 if ct == 'tp' and is_ins and is_ins_truth else tp_ins
                tp_del = tp_del + 1 if ct == 'tp' and is_del and is_del_truth else tp_del

                # fp_snp_truth = fp_snp_truth + 1 if ct == 'tp' and is_snp_truth and not is_snp else fp_snp_truth
                # fp_ins_truth = fp_ins_truth + 1 if ct == 'tp' and is_ins_truth and not is_ins else fp_ins_truth
                # fp_del_truth = fp_del_truth + 1 if ct == 'tp' and is_del_truth and not is_del else fp_del_truth

                fp_snp = fp_snp + 1 if is_snp and ct == 'fp' else fp_snp
                fp_ins = fp_ins + 1 if is_ins and ct == 'fp' else fp_ins
                fp_del = fp_del + 1 if is_del and ct == 'fp' else fp_del

                fn_snp = fn_snp + 1 if is_snp and ct == 'fn' else fn_snp
                fn_ins = fn_ins + 1 if is_ins and ct == 'fn' else fn_ins
                fn_del = fn_del + 1 if is_del and ct == 'fn' else fn_del

    print(''.join([x.ljust(15) for x in ['tp_snp', 'tp_ins', 'tp_del', 'fp_snp', 'fp_ins', 'fp_del', 'fp_snp_truth', 'fp_ins_truth', 'fp_del_truth', 'fn_snp', 'fn_ins', 'fn_del']]))
    print(''.join([str(x).ljust(15) for x in [tp_snp, tp_ins, tp_del, fp_snp, fp_ins, fp_del,fp_snp_truth,fp_ins_truth, fp_del_truth, fn_snp, fn_ins, fn_del]]))

    # cal f1-score
    snp_pre, snp_rec, snp_f1 = cal_metrics(tp=tp_snp, fp=fp_snp+fp_snp_truth, fn=fn_snp)
    ins_pre, ins_rec, ins_f1 = cal_metrics(tp=tp_ins, fp=fp_ins+fp_ins_truth, fn=fn_ins)
    del_pre, del_rec, del_f1 = cal_metrics(tp=tp_del, fp=fp_del+fp_del_truth, fn=fn_del)

    print(''.join([x.ljust(15) for x in ['Variant type', 'Precision', 'Recall', 'F1-score']]))
    print(''.join([str(x).ljust(15) for x in ['SNP', snp_pre, snp_rec, snp_f1]]))
    print(''.join([str(x).ljust(15) for x in ['INS', ins_pre, ins_rec, ins_f1]]))
    print(''.join([str(x).ljust(15) for x in ['DEL', del_pre, del_rec, del_f1]]))


        # if output_logits:
        #     prediction_pro, normal_logit, tumor_logit = prediction
        #     normal_logits = np.concatenate([normal_logits, normal_logit], axis=0)
        #     tumor_logits = np.concatenate([tumor_logits, tumor_logit], axis=0)
        # labels = np.concatenate([labels, val_label[0]], axis = 0)
        # if val_total > 100000:
        #     break

    # if output_dir and os.path.exists(output_dir):
    #     if output_logits:
    #         np.save(os.path.join(output_dir, 'normal'), normal_logits)
    #         np.save(os.path.join(output_dir, 'tumor'), tumor_logits)
    #     np.save(os.path.join(output_dir, 'label'), labels)
    #     np.save(os.path.join(output_dir, 'prediction'), predictions)
    logging.info("[INFO] Total val data size/{}".format(val_total))

    fp_vcf.close()
    fn_vcf.close()
    tp_vcf.close()
    # validate_dataset = validate_dataset if add_validation_dataset else None
    # if args.chkpnt_fn is not None:
    #     model.load_weights(args.chkpnt_fn)
    #
    # train_history = model.fit(x=train_dataset,
    #                           epochs=max_epoch,
    #                           validation_data=validate_dataset,
    #                           callbacks=[early_stop_callback, model_save_callbakck],
    #                           verbose=1,
    #                           shuffle=False)
    #
    # for table_dataset in table_dataset_list:
    #     table_dataset.close()
    #
    # for table_dataset in validate_table_dataset_list:
    #     table_dataset.close()
    #
    # # show the parameter set with the smallest validation loss
    # if 'val_loss' in train_history.history:
    #     best_validation_epoch = np.argmin(np.array(train_history.history["val_loss"])) + 1
    #     logging.info("[INFO] Best validation loss at epoch: %d" % best_validation_epoch)
    # else:
    #     best_train_epoch = np.argmin(np.array(train_history.history["loss"])) + 1
    #     logging.info("[INFO] Best train loss at epoch: %d" % best_train_epoch)


def main():
    parser = ArgumentParser(description="Train a Clair3 model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bin_fn', type=str, default="", required=True,
                        help="Binary tensor input generated by Tensor2Bin.py, support multiple bin readers using pytables")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a model to resume training or for fine-tuning")

    parser.add_argument('--ochk_prefix', type=str, default=None,
                        help="Prefix for model output after each epoch")

    # options for advanced users
    parser.add_argument('--maxEpoch', type=int, default=None,
                        help="Maximum number of training epochs")

    parser.add_argument('--learning_rate', type=float, default=None,
                        help="Set the initial learning rate, default: %(default)s")


    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="Define training samples to be excluded")

    parser.add_argument('--unified_vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--output_logits', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Prefix for model output after each epoch")

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--use_resnet', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--add_writer', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--use_siam', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--add_contrastive', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--test_all_pos', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--predict_germline', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    vgrp = parser.add_mutually_exclusive_group()
    vgrp.add_argument('--random_validation', action='store_true',
                        help="Use random sample of dataset for validation, default: %(default)s")

    vgrp.add_argument('--validation_fn', type=str, default=None,
                        help="Binary tensor input for use in validation: %(default)s")


    args = parser.parse_args()


    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    predict_model(args)


if __name__ == "__main__":
    main()
# /autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/python3 /mnt/bal36/zxzheng/somatic/Clair-somatic/clair-somatic.py Predict --bin_fn /mnt/bal36/zxzheng/somatic/hifi/all_add_hete/build/bins --ochk_prefix /mnt/bal36/zxzheng/somatic/hifi/all_add_hete/train/test --platform ont