import logging
import random
import numpy as np
from argparse import ArgumentParser, SUPPRESS
# import tensorflow_addons as tfa
# import tensorflow as tf
import tables
import os
import sys
from itertools import accumulate
from time import time
# import clair_somatic.model as model_path
from shared.utils import str2bool
import shared.param as param
logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '64'
os.environ['NUMEXPR_NUM_THREADS'] = '8'
from tqdm import tqdm
import os
import random
import zipfile

# import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
# from linformer import Linformer
# from PIL import Image
# from sklearn.model_selection import train_test_split
from torch.optim.lr_scheduler import StepLR
from torch.utils.data import DataLoader, Dataset
# from torchvision import datasets, transforms
# from tqdm.notebook import tqdm

from torchsummary import summary
from clair_somatic.cvt import CvT

# reuqired package  torchsummary, tqdm tables,  einops

batch_size = 64
epochs = 20
lr = 3e-5
gamma = 0.7
seed = 42
device = 'cuda'

random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True


# import glob
# from itertools import chain
# import os
# import random
# import zipfile
#
# import matplotlib.pyplot as plt
# import numpy as np
# # import pandas as pd
# import torch
# import torch.nn as nn
# import torch.nn.functional as F
# import torch.optim as optim
# from linformer import Linformer
# from PIL import Image
# from sklearn.model_selection import train_test_split
# from torch.optim.lr_scheduler import StepLR
# from torch.utils.data import DataLoader, Dataset
# from torchvision import datasets, transforms
# from tqdm.notebook import tqdm
#
# from vit_pytorch.efficient import ViT
#
# model = ViT(
#     dim=128,
#     image_size=224,
#     patch_size=32,
#     num_classes=2,
#     transformer=efficient_transformer,
#     channels=3,
# ).to(device)
#


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

#
# class FocalLoss(tf.keras.losses.Loss):
#     """
#     updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
#     more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
#     """
#
#     def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
#         super(FocalLoss, self).__init__()
#         self.gamma = gamma
#         self.cls_weights = None
#         if effective_label_num is not None:
#             task_label_num = get_label_task(effective_label_num, label_shape_cum, task)
#             cls_weights = cal_class_weight(task_label_num, len(task_label_num))
#             cls_weights = tf.constant(cls_weights, dtype=tf.float32)
#             cls_weights = tf.expand_dims(cls_weights, axis=0)
#             self.cls_weights = cls_weights
#
#     def call(self, y_true, y_pred):
#         y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
#         cross_entropy = -y_true * tf.math.log(y_pred)
#         weight = ((1 - y_pred) ** self.gamma) * y_true
#         FCLoss = cross_entropy * weight
#         if self.cls_weights is not None:
#             FCLoss = FCLoss * self.cls_weights
#         reduce_fl = tf.reduce_sum(FCLoss, axis=-1)
#         return reduce_fl
#
# class BinaryCrossentropy(tf.keras.losses.Loss):
#     """
#     updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
#     more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
#     """
#
#     def __init__(self):
#         super(BinaryCrossentropy, self).__init__()
#
#     def call(self, y_true, y_pred):
#         sigmoids = tf.nn.sigmoid_cross_entropy_with_logits(labels=y_true, logits=y_pred)
#         sigmoids_loss = tf.reduce_mean(sigmoids)
#         return sigmoids_loss


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
            # print ("Sum:{}".format(np.sum(np.array(all_shuffle_chunk_list[:offset_idx]))))
            return np.array(all_shuffle_chunk_list[:offset_idx]), np.array(all_shuffle_chunk_list[offset_idx + 1:])
        else:
            total_size += chunk_num * chunk_size
            offset_idx += chunk_num


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

# def compute_euclidean_distance(x, y):
#     """
#     Computes the euclidean distance between two tensorflow variables
#     """
#
#     d = tf.reduce_sum(tf.square(tf.sub(x, y)),1)
#     return d


#
# class ContrastiveLoss(tf.keras.losses.Loss):
#
#     """
#     Compute the contrastive loss as in
#
#
#     L = 0.5 * Y * D^2 + 0.5 * (Y-1) * {max(0, margin - D)}^2
#     Y=1: same class(similar), Y=0: different class
#     **Parameters**
#      left_feature: First element of the pair
#      right_feature: Second element of the pair
#      label: Label of the pair (0 or 1)
#      margin: Contrastive margin
#
#     **Returns**
#      Return the loss operation
#
#     """
#
#     def __init__(self, margin=1):
#         super(ContrastiveLoss, self).__init__()
#         self.margin = margin
#
#     # def call(self, y_true, y_pred):
#     #     label = tf.argmax(y_true, axis=1)
#     #     label = tf.cast(label, tf.float32)
#     #
#     #     d_sqrt = tf.sqrt(y_pred)
#     #     first_part = tf.multiply(1.0 - label, y_pred)  # (Y-1)*(d)
#     #
#     #     max_part = tf.square(tf.maximum(self.margin - d_sqrt, 0))
#     #     second_part = tf.multiply(label, max_part)  # (Y) * max(margin - d, 0)
#     #
#     #     loss = 0.5 * tf.reduce_sum(first_part + second_part, axis=-1)
#     #     return loss
#     def call(self, y_true, y_pred):
#         label = tf.argmax(y_true, axis=1)
#         label = tf.cast(label, tf.float32)
#
#         d_sqrt = tf.sqrt(y_pred)
#         first_part = tf.multiply(1.0 - label, y_pred)  # (Y-1)*(d)
#
#         max_part = tf.square(tf.maximum(self.margin - d_sqrt, 0))
#         second_part = tf.multiply(label, max_part)  # (Y) * max(margin - d, 0)
#
#         loss = 0.5 * tf.reduce_sum(first_part + second_part, axis=-1)
#         return loss

def train_model(args):
    apply_focal_loss = False

    platform = args.platform
    ctg_name_string = args.ctgName
    chkpnt_fn = args.chkpnt_fn
    ochk_prefix = args.ochk_prefix
    add_validation_dataset = args.random_validation or (args.validation_fn is not None)
    validation_fn = args.validation_fn
    ctg_name_list = ctg_name_string.split(',')  if ctg_name_string is not None else []
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()

    if ochk_prefix and not os.path.exists(ochk_prefix):
        from subprocess import run
        output = run('mkdir {}'.format(ochk_prefix), shell=True)
        print("[INFO] Model path empty, mkdir folder")
    model = CvT(
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
        apply_softmax = apply_focal_loss
    ).to(device)

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
    bin_list = bin_list[:3]
    # bin_list = bin_list[:10]
    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))
    # total_data_size = 0
    # table_dataset_list = []
    # validate_table_dataset_list = []
    # chunk_offset = np.zeros(len(bin_list), dtype=int)

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
    # for bin_idx, bin_file in enumerate(bin_list):
    #     table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
    #     validate_table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
    #     table_dataset_list.append(table_dataset)
    #     validate_table_dataset_list.append(validate_table_dataset)
    #     chunk_num = (len(table_dataset.root.label) - batch_size) // chunk_size
    #     data_size = int(chunk_num * chunk_size)
    #     chunk_offset[bin_idx] = chunk_num
    #     total_data_size += data_size

    # train_data_size = total_data_size * param.trainingDatasetPercentage
    # train_data_size = int(train_data_size // chunk_size) * chunk_size
    # validate_data_size = int((total_data_size - train_data_size) // chunk_size) * chunk_size
    # train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_data_size, chunk_size)

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
    def DataGenerator(x, shuffle_chunk_list, train_flag=True):

        batch_num = len(shuffle_chunk_list) // chunks_per_batch
        # normal_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        # tumor_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        input_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        label = np.empty((batch_size, param.label_size), np.float32)

        for epoch in range(epochs):
            random_start_position = np.random.randint(0, batch_size) if train_flag else 0
            if train_flag:
                np.random.shuffle(shuffle_chunk_list)
            for batch_idx in range(batch_num):

                for chunk_idx in range(chunks_per_batch):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunks_per_batch + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id
                    # normal_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.normal_matrix[
                    #         random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    #
                    # tumor_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.tumor_matrix[
                    #         random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    #
                    input_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.input_matrix[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]

                # input_matrix = np.concatenate((normal_matrix, tumor_matrix), axis=1)
                # if add_contrastive:
                label_for_normal = [[0,1] if np.argmax(item) == 1 else [1,0] for item in label]
                label_for_tumor = [[0,0,1] if np.argmax(item) == 2 else ([0, 1,0] if np.argmax(item) == 1 else [1,0,0]) for item in label]
                label_for_normal = np.array(label_for_normal, dtype=np.float32)
                label_for_tumor = np.array(label_for_tumor, dtype=np.float32)

                input_tensor = torch.from_numpy(np.transpose(input_matrix, (0,3,1,2))).to('cuda')
                # tumor_tensor = torch.from_numpy(tumor_matrix)
                label_tensor = torch.from_numpy(label_for_tumor).to('cuda')
                yield input_tensor, label_tensor

    print (summary(model, input_size=(param.channel_size,param.max_depth,param.no_of_positions)))
    train_dataset_loder = DataGenerator(table_dataset_list, train_shuffle_chunk_list, True)
    validate_dataset_loder = DataGenerator(validate_table_dataset_list, validate_shuffle_chunk_list, False)

    criterion = nn.CrossEntropyLoss()
    # criterion = FocalLoss()
    # optimizer
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    # scheduler
    scheduler = StepLR(optimizer, step_size=1, gamma=gamma)
    
    train_steps = train_data_size // batch_size
    validate_steps = validate_data_size // batch_size
    print("[INFO] The size of dataset: {}".format(train_data_size))
    print("[INFO] The training batch size: {}".format(batch_size))
    print("[INFO] The training learning_rate: {}".format(learning_rate))
    print("[INFO] The output model folder: {}".format(ochk_prefix))

    print('[INFO] Train steps:{}'.format(train_steps))
    print('[INFO] Validate steps:{}'.format(validate_steps))
    for epoch in range(epochs):
        epoch_loss = 0
        start_time = time()
        tumor_fp, tumor_tp, tumor_fn = 0,0,0
        val_tumor_fp, val_tumor_tp, val_tumor_fn = 0,0,0
        ref_fp, ref_tp, ref_fn = 0,0,0

        t = tqdm(enumerate(train_dataset_loder), total=train_steps, position=0, leave=True)
        for batch_idx, (data, label) in t:
            t.set_description('EPOCH {}'.format(epoch))
            data = data.to(device)
            label = label.to(device)
            current_batch_size = len(label)

            output_logit = model(data.float())
            y_truth = torch.argmax(label, axis=1)
            if apply_focal_loss:
                # y_pred, y_true
                loss = criterion(input=output_logit, target=label)
            else:
                # CE
                loss = criterion(output_logit, y_truth)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            y_truth = y_truth.cpu().numpy()
            y_pred = output.argmax(dim=1).cpu().numpy()
            arg_index = 2
            tumor_fp += sum([True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)])
            tumor_fn += sum([True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)])
            tumor_tp += sum([True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)])

            if batch_idx + 1 == train_steps:
                break

            epoch_loss += loss
            el = epoch_loss.detach().cpu().numpy()
            t.set_postfix({'loss': el, 'tp': tumor_tp, 'fp': tumor_fp, 'fn': tumor_fn,})
            t.update(1)


        # with torch.no_grad():
        #     epoch_val_loss = 0
        #     for data, label in validate_dataset_loder:
        #         data = data.to(device)
        #         label = label.to(device)
        #         current_batch_size = len(label)
        #
        #         val_output = model(data.float())
        #         val_loss = criterion(val_output, label)
        #
        #         acc = (val_output.argmax(dim=1) == label).float().mean()
        #         epoch_val_accuracy += acc / current_batch_size
        #         epoch_val_loss += val_loss / current_batch_size
        #
        #     print(f"Epoch : {epoch+1} - loss : {epoch_loss:.4f} - acc: {epoch_accuracy:.4f} - val_loss : {epoch_val_loss:.4f} - val_acc: {epoch_val_accuracy:.4f}\n")

        save_path = os.path.join(ochk_prefix, "{}.pkl".format(epoch)) if ochk_prefix is not None else "{}.pkl".format(epoch)
        torch.save(model, save_path)
        # model = torch.load("{}.pkl".format(epoch))


    for table_dataset in table_dataset_list:
        table_dataset.close()

    for table_dataset in validate_table_dataset_list:
        table_dataset.close()

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

    parser.add_argument('--learning_rate', type=float, default=1e-3,
                        help="Set the initial learning rate, default: %(default)s")


    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")

    parser.add_argument('--use_siam', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--add_contrastive', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--ctgName', type=str, default=None,
                        help="Define training samples to be excluded")
    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    # mutually-incompatible validation options
    vgrp = parser.add_mutually_exclusive_group()
    vgrp.add_argument('--random_validation', action='store_true',
                        help="Use random sample of dataset for validation, default: %(default)s")

    vgrp.add_argument('--validation_fn', type=str, default=None,
                        help="Binary tensor input for use in validation: %(default)s")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    train_model(args)


if __name__ == "__main__":
    main()
