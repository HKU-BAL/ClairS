import logging
import random
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import tensorflow_addons as tfa
import tensorflow as tf
import tables
import os
import sys
from itertools import accumulate
from collections import defaultdict

import clair_somatic.model as model_path
from shared.utils import str2bool
from shared.vcf import VcfReader

logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '64'
os.environ['NUMEXPR_NUM_THREADS'] = '8'


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


class FocalLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.cls_weights = None
        if effective_label_num is not None:
            task_label_num = get_label_task(effective_label_num, label_shape_cum, task)
            cls_weights = cal_class_weight(task_label_num, len(task_label_num))
            cls_weights = tf.constant(cls_weights, dtype=tf.float32)
            cls_weights = tf.expand_dims(cls_weights, axis=0)
            self.cls_weights = cls_weights

    def call(self, y_true, y_pred):
        y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        cross_entropy = -y_true * tf.math.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight
        if self.cls_weights is not None:
            FCLoss = FCLoss * self.cls_weights
        reduce_fl = tf.reduce_sum(FCLoss, axis=-1)
        return reduce_fl

class BinaryCrossentropy(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self):
        super(BinaryCrossentropy, self).__init__()

    def call(self, y_true, y_pred):
        sigmoids = tf.nn.sigmoid_cross_entropy_with_logits(labels=y_true, logits=y_pred)
        sigmoids_loss = tf.reduce_mean(sigmoids)
        return sigmoids_loss


def get_chunk_list(chunk_offset, train_data_size, chunk_size):
    """
    get chunk list for training and validation data. we will randomly split training and validation dataset,
    all training data is directly acquired from various tensor bin files.

    """
    all_shuffle_chunk_list = []
    total_size = 0
    offset_idx = 0
    for bin_idx, chunk_num in enumerate(chunk_offset):
        all_shuffle_chunk_list += [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
    np.random.seed(0)
    np.random.shuffle(all_shuffle_chunk_list)  # keep the same random validate dataset
    for bin_idx, chunk_num in enumerate(chunk_offset):
        if chunk_num * chunk_size + total_size >= train_data_size:
            chunk_num = (train_data_size - total_size) // chunk_size
            offset_idx += chunk_num
            print ("Sum:{}".format(np.sum(np.array(all_shuffle_chunk_list[:offset_idx]))))
            return np.array(all_shuffle_chunk_list[:offset_idx]), np.array(all_shuffle_chunk_list[offset_idx + 1:])
        else:
            total_size += chunk_num * chunk_size
            offset_idx += chunk_num


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
    depth, alt_str, af_str = alt_info.split('-')
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

def predict_model(args):
    platform = args.platform
    pileup = args.pileup
    use_siam = args.use_siam
    add_contrastive = args.add_contrastive
    add_indel_length = args.add_indel_length
    output_logits = args.output_logits
    output_dir = args.output_dir
    ctg_name_string = args.ctgName
    predict_germline = args.predict_germline
    ctg_name_list = ctg_name_string.split(',')  if ctg_name_string is not None else []
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    add_validation_dataset = True
    test_all_pos = args.test_all_pos
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""
    is_gpu_available = tf.test.is_gpu_available()
    use_tf_dataset = False
    unified_vcf_fn = args.unified_vcf_fn

    unified_vcf_reader_dict = defaultdict()
    for ctg_name in ctg_name_list:
        unified_vcf_reader = VcfReader(vcf_fn=unified_vcf_fn + ctg_name, ctg_name=ctg_name, is_var_format=False)
        unified_vcf_reader.read_vcf()
        unified_vcf_reader_dict[ctg_name] = unified_vcf_reader.variant_dict

    if not is_gpu_available:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
    else:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        tf.config.experimental.set_virtual_device_configuration(gpus[0], [
            tf.config.experimental.VirtualDeviceConfiguration(memory_limit=2048)])
    if pileup:
        import shared.param_p as param
        model = model_path.Clair3_P()
    else:
        import shared.param as param
        # if use_siam:
        #     model = model_path.Clair3_Siam(add_indel_length=add_indel_length, add_contrastive=add_contrastive, output_logits=output_logits)
        # else:
        model = model_path.Clair3_F(add_indel_length=add_indel_length)

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    label_size, label_shape = param.label_size, param.label_shape
    label_shape_cum = list(accumulate(label_shape))
    batch_size, chunk_size = param.predictBatchSize, param.test_chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    task_num = 1
    # TensorShape = (tuple(tf.TensorShape([None] + tensor_shape) for _ in range(2)),
    #     tuple(tf.TensorShape([None, label_shape[task]]) for task in range(task_num)))
    #
    # TensorDtype = ((tf.int32,tf.int32), tuple(tf.float32 for _ in range(task_num)))

    Val_TensorShape = (tuple(tf.TensorShape([None] + tensor_shape) for _ in range(2)),
        tuple(tf.TensorShape([None, label_shape[task]]) for task in range(task_num)), tf.TensorShape([None] + [1]), tf.TensorShape([None] + [1]), tf.TensorShape([None] + [1]))

    Val_TensorDtype = ((tf.int32,tf.int32), tuple(tf.float32 for _ in range(task_num)), tf.string, tf.string, tf.string)

    bin_list = os.listdir(args.bin_fn)
    # default we exclude sample hg003 and all chr20 for training
    bin_list = [f for f in bin_list if pass_chr(f, ctg_name_list) and not exist_file_prefix(exclude_training_samples, f)]
    # bin_list = [f for f in bin_list if f.endswith('tp')]
    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))
    total_data_size = 0
    table_dataset_list = []
    validate_table_dataset_list = []
    chunk_offset = np.zeros(len(bin_list), dtype=int)

    for bin_idx, bin_file in enumerate(bin_list):
        table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        validate_table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        table_dataset_list.append(table_dataset)
        validate_table_dataset_list.append(validate_table_dataset)
        chunk_num = (len(table_dataset.root.label) - batch_size) // chunk_size
        data_size = int(chunk_num * chunk_size)
        chunk_offset[bin_idx] = chunk_num
        total_data_size += data_size
    if test_all_pos:
        train_data_size = 0
    else:
        train_data_size = total_data_size * param.trainingDatasetPercentage
    validate_data_size = int((total_data_size - train_data_size) // chunk_size) * chunk_size
    train_data_size = int(train_data_size // chunk_size) * chunk_size
    train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_data_size, chunk_size)

    def DataGenerator(x, data_size, shuffle_chunk_list, train_flag=True):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """

        chunk_iters = batch_size // chunk_size
        batch_num = data_size // batch_size
        normal_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        tumor_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        label = np.empty((batch_size, param.label_size), np.float32)
        positions = np.empty((batch_size, 1), 'S2000')
        alt_infos = np.empty((batch_size, 1), 'S2000')
        ctg_name_infos = np.empty((batch_size, 1), 'S20')

        random_start_position = np.random.randint(0, batch_size) if train_flag else 0
        if train_flag:
            np.random.shuffle(shuffle_chunk_list)
        for batch_idx in range(batch_num):
            for chunk_idx in range(chunk_iters):
                offset_chunk_id = shuffle_chunk_list[batch_idx * chunk_iters + chunk_idx]
                bin_id, chunk_id = offset_chunk_id
                normal_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.normal_matrix[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

                tumor_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.tumor_matrix[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

                label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

                if not train_flag:
                    positions[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    alt_infos[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.alt_info[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    contig_infos = x[bin_id].filename.split('/')[-1].split('_')[1].split('.')[0]
                    ctg_name_infos[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = [[contig_infos]] * chunk_size

            if train_flag:
                yield (normal_matrix, tumor_matrix), (label[:,:label_shape_cum[0]],)
            else:
                yield (normal_matrix, tumor_matrix), (label[:, :label_shape_cum[0]],), positions, alt_infos, ctg_name_infos
            # else:
            #     yield (normal_matrix, tumor_matrix), (label[:,:label_shape_cum[0]])

    def DataGenerator_all(x, data_size, shuffle_chunk_list, train_flag=True):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """

        for dataset in x:
            data_size = len(dataset.root.label)
            batch_num = data_size // batch_size if data_size % batch_size == 0 else data_size // batch_size + 1

            for batch_idx in range(batch_num):
                current_batch_size = batch_size if batch_idx != batch_num -1 else data_size % batch_size
                normal_matrix = np.empty([current_batch_size] + tensor_shape, np.int32)
                tumor_matrix = np.empty([current_batch_size] + tensor_shape, np.int32)
                label = np.empty((current_batch_size, param.label_size), np.float32)
                positions = np.empty((current_batch_size, 1), 'S2000')
                alt_infos = np.empty((current_batch_size, 1), 'S2000')
                ctg_name_infos = np.empty((current_batch_size, 1), 'S20')

                normal_matrix = np.array(dataset.root.normal_matrix[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                tumor_matrix = np.array(dataset.root.tumor_matrix[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                label = np.array(dataset.root.label[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                positions = np.array(dataset.root.position[batch_idx * batch_size:(batch_idx + 1) * batch_size])
                alt_infos = np.array(dataset.root.alt_info[batch_idx * batch_size:(batch_idx + 1) * batch_size])

                contig_infos = dataset.filename.split('/')[-1].split('_')[1].split('.')[0]
                ctg_name_infos = np.array([[contig_infos]] * current_batch_size)

                yield (normal_matrix, tumor_matrix), (label[:, :label_shape_cum[0]]), positions, alt_infos, ctg_name_infos

    # train_dataset = tf.data.Dataset.from_generator(
    #     lambda: DataGenerator(table_dataset_list, train_data_size, train_shuffle_chunk_list, True), TensorDtype,
    #     TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    if test_all_pos:
        validate_dataset = DataGenerator_all(validate_table_dataset_list, validate_data_size, validate_shuffle_chunk_list, False)
    else:
        validate_dataset = tf.data.Dataset.from_generator(
            lambda: DataGenerator(validate_table_dataset_list, validate_data_size, validate_shuffle_chunk_list, False),
            Val_TensorDtype,
            Val_TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    total_steps = max_epoch * train_data_size // batch_size

    #RectifiedAdam with warmup start
    optimizer = tfa.optimizers.Lookahead(tfa.optimizers.RectifiedAdam(
        lr=learning_rate,
        total_steps=total_steps,
        warmup_proportion=0.1,
        min_lr=learning_rate*0.75,
    ))
    # optimizer = tf.optimizers.Adam()
    # loss_func = [FocalLoss(label_shape_cum, task, effective_label_num) for task in range(task_num)]
    loss_func = [BinaryCrossentropy() for task in range(task_num)]
    loss_task = {"output_{}".format(task + 1): loss_func[task] for task in range(task_num)}
    metrics = {"output_{}".format(task + 1): tfa.metrics.F1Score(num_classes=label_shape[task], average='micro') for
               task in range(task_num)}

    model.compile(
        loss=loss_task,
        metrics=metrics,
        optimizer=optimizer
    )
    early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, mode="min")
    model_save_callbakck = tf.keras.callbacks.ModelCheckpoint(os.path.join(ochk_prefix, "{epoch}") if ochk_prefix else "{epoch}", period=1, save_weights_only=False)

    # Use first 20 element to initialize tensorflow model using graph mode
    output = model((np.array(table_dataset_list[0].root.normal_matrix[:20]),np.array(table_dataset_list[0].root.tumor_matrix[:20])))
    logging.info(model.summary(print_fn=logging.info))

    if args.chkpnt_fn is not None:
        model.load_weights(args.chkpnt_fn)

    val_total = 0
    predictions = np.empty(shape=(0, label_size))
    normal_logits = np.empty(shape=(0, 64))
    tumor_logits = np.empty(shape=(0, 64))
    labels = np.empty(shape=(0, label_size))

    prefix = 'germline_' if predict_germline else ""
    fp_vcf = open(os.path.join(output_dir, prefix +'fp.vcf'), 'w')
    fn_vcf = open(os.path.join(output_dir, prefix+'fn.vcf'), 'w')
    tp_vcf = open(os.path.join(output_dir, prefix+'tp.vcf'), 'w')

    tp_snp, tp_ins, tp_del, fp_snp, fp_ins, fp_del, fn_snp, fn_ins, fn_del, fp_snp_truth, fp_ins_truth, fp_del_truth = 0,0,0,0,0,0,0,0,0,0,0,0
    tp_count = 0
    for val_data, val_label, val_positions, val_alt_infos, val_ctg_name_infos in validate_dataset:
        prediction = model.predict_on_batch(val_data)
        # prediction = tf.nn.sigmoid(prediction[0]).numpy()
        try:
            prediction = prediction[0].numpy()
        except:
            prediction = prediction[0]
        val_total += len(prediction)
        # predictions = np.concatenate([predictions, prediction], axis = 0)
        y_pred = np.argmax(prediction, axis=1)
        y_truth = np.argmax(val_label, axis=1)
        #np where to get index
        arg_index = 1 if predict_germline else 2
        fp = [True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)]
        fn = [True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)]
        tp = [True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)]
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

        fp_pos_array = [(item[0],item[-1]) for item in fp_pos_array]
        fn_pos_array = [(item[0],item[-1]) for item in fn_pos_array]
        tp_pos_array = [(item[0],item[-1]) for item in tp_pos_array]

        all_pos_array = [fp_pos_array, fn_pos_array, tp_pos_array]
        all_chr_array = [fp_chr_array, fn_chr_array, tp_chr_array]
        all_alt_array = [fp_alt_array, fn_alt_array, tp_alt_array]
        all_vcf_fn = [fp_vcf, fn_vcf, tp_vcf]
        candidate_types = ['fp', 'fn', 'tp']
        for pos_array, chr_array, vcf_fn, alt_array, ct in zip(all_pos_array, all_chr_array, all_vcf_fn, all_alt_array, candidate_types):
            # print(len(pos_array), len(chr_array), len(alt_array), len(ct))
            for (pos, variant_type), contig, alt_info in zip(pos_array, chr_array, alt_array):
                vcf_format = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF:VT\t%s:%d:%d:%.4f:%s" % (
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
                    0.5,
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

    parser.add_argument('--validation_dataset', action='store_true',
                        help="Use validation dataset when training, default: %(default)s")

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

    parser.add_argument('--use_siam', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--add_contrastive', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--test_all_pos', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--predict_germline', action='store_true',
                        help=SUPPRESS)
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    predict_model(args)


if __name__ == "__main__":
    main()
# /autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/python3 /mnt/bal36/zxzheng/somatic/Clair-somatic/clair-somatic.py Predict --bin_fn /mnt/bal36/zxzheng/somatic/hifi/all_add_hete/build/bins --ochk_prefix /mnt/bal36/zxzheng/somatic/hifi/all_add_hete/train/test --platform ont