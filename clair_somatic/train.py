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
import torch.nn.functional as F
from torch.autograd import Variable

from tqdm import tqdm

from shared.utils import str2bool
import shared.param as param
import clair_somatic.model as model_path

# reuqired package  torchsummary, tqdm tables,  einops
logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '256'
os.environ['NUMEXPR_NUM_THREADS'] = '32'
gamma = 0.7
seed = 100
random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
from torch.autograd import Variable
#
# class FocalLoss(nn.Module):
#     def __init__(self, gamma=2, alpha=None, size_average=True):
#         super(FocalLoss, self).__init__()
#         self.gamma = gamma
#         self.alpha = alpha
#         if isinstance(alpha,(float,int)): self.alpha = torch.Tensor([alpha,1-alpha])
#         if isinstance(alpha,list): self.alpha = torch.Tensor(alpha)
#         self.size_average = size_average
#
#     def forward(self, input, target):
#         if input.dim() > 2:
#             input = input.view(input.size(0),input.size(1),-1)  # N,C,H,W => N,C,H*W
#             input = input.transpose(1,2)    # N,C,H*W => N,H*W,C
#             input = input.contiguous().view(-1,input.size(2))   # N,H*W,C => N*H*W,C
#         target = torch.argmax(target, axis=1)
#         target = target.view(-1,1)
#
#         logpt = F.log_softmax(input,dim=-1)
#         logpt = logpt.gather(1,target)
#         logpt = logpt.view(-1)
#         pt = Variable(logpt.data.exp())
#
#         if self.alpha is not None:
#             if self.alpha.type()!=input.data.type():
#                 self.alpha = self.alpha.type_as(input.data)
#             at = self.alpha.gather(0,target.data.view(-1))
#             logpt = logpt * Variable(at)
#
#         loss = -1 * (1-pt)**self.gamma * logpt
#         if self.size_average:
#             return loss.mean()
#         else:
#             return loss.sum()

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
        y_pred = torch.nn.functional.softmax(y_pred, dim=1)
        y_pred = torch.clamp(y_pred, min=1e-9, max=1-1e-9)
        cross_entropy = -y_true * torch.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight

        reduce_fl = torch.mean(torch.sum(FCLoss,dim=1))
        return reduce_fl


class AFLoss(nn.Module):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, alpha=None, gamma=2):
        super(AFLoss, self).__init__()
        self.gamma = gamma

    def forward(self, input, target):
        y_pred, y_true = input, target
        y_pred = torch.nn.functional.softmax(y_pred, dim=1)
        y_pred = torch.clamp(y_pred, min=1e-9, max=1-1e-9)

        cross_entropy = -y_true * torch.log(y_pred)
        reduce_fl = torch.mean(torch.sum(cross_entropy, dim=1))
        # cross_entropy = -y_true * torch.log(y_pred)
        # weight = ((1 - y_pred) ** self.gamma) * y_true
        # FCLoss = cross_entropy * weight
        # reduce_fl = torch.mean(torch.sum(FCLoss,dim=1))

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
        if ctg_name + '.' in fn:
            return True
    return False

def train_model(args):
    apply_focal_loss = param.apply_focal_loss
    discard_germline = param.discard_germline
    add_l2_regulation_loss = param.add_l2_regulation_loss
    l2_regularization_lambda = param.l2_regularization_lambda

    use_resnet = args.use_resnet
    platform = args.platform
    ctg_name_string = args.ctg_name
    chkpnt_fn = args.chkpnt_fn
    ochk_prefix = args.ochk_prefix
    add_writer = args.add_writer

    add_validation_dataset = args.random_validation or (args.validation_fn is not None)
    validation_fn = args.validation_fn
    ctg_name_list = ctg_name_string.split(',') if ctg_name_string is not None else []
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()

    if ochk_prefix and not os.path.exists(ochk_prefix):
        output = run('mkdir -p {}'.format(ochk_prefix), shell=True)
        print("[INFO] Model path empty, create folder:{}".format(ochk_prefix))

    if add_writer:
        writer = SummaryWriter("{}/log".format(ochk_prefix))
    device = 'cpu'
    if torch.cuda.is_available():
        # gpu_id = torch.cuda.current_device()
        # total_memory = torch.cuda.get_device_properties(gpu_id).total_memory
        # reserved = torch.cuda.memory_reserved(gpu_id)
        # allocated = torch.cuda.memory_allocated(gpu_id)
        # free_memory = (reserved - allocated) / (1024 ** 3)  # free inside reserved
        # if free_memory >= 1:
        device = 'cuda'
    apply_softmax = False if apply_focal_loss else True
    if args.pileup:
        model = model_path.bigru(apply_softmax=apply_softmax, num_classes=2 if discard_germline else 3).to(device)
        channel_size = param.pileup_channel_size
        pileup_tensor_shape = [param.no_of_positions, channel_size * 2] # normal and tumor
        tensor_shape = pileup_tensor_shape
        input = torch.ones(size=[100] + tensor_shape).to(device)

    else:
        #use torchinfo
        model = model_path.CvT(
            num_classes=2 if discard_germline else 3,
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
            apply_softmax=apply_softmax
        ).to(device)
        input = torch.ones(size=(100, param.channel_size, param.max_depth, param.no_of_positions)).to(device)
        tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape

        if use_resnet:
            model = model_path.ResNet().to(device)

    if chkpnt_fn is not None:
        model = torch.load(chkpnt_fn)

    output = model(input)

    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    # if args.pileup:
    #     batch_size = 200
    assert batch_size % chunk_size == 0
    chunks_per_batch = batch_size // chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    bin_list = os.listdir(args.bin_fn)

    bin_list = [f for f in bin_list if pass_chr(f, ctg_name_list) and not exist_file_prefix(exclude_training_samples, f)]
    failed_bin_set = set()
    for bin_file in bin_list:
        try:
            table = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
            table.close()
        except:
            print("[WARNING] {} cannot open!".format(bin_file))
            failed_bin_set.add(bin_file)
    bin_list = [f for f in bin_list if f not in failed_bin_set]
    if len(bin_list) == 0:
        print("[ERROR] Cannot find ant binary for model training")
        return

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


    def DataGenerator(x, shuffle_chunk_list, train_flag=True):

        batch_num = len(shuffle_chunk_list) // chunks_per_batch
        input_matrix = np.empty([batch_size] + tensor_shape, np.float32)
        label = np.empty((batch_size, param.label_size), np.float32)

        for epoch in range(max_epoch):
            random_start_position = np.random.randint(0, batch_size) if train_flag else 0
            if train_flag:
                np.random.shuffle(shuffle_chunk_list)
            for batch_idx in range(batch_num):

                for chunk_idx in range(chunks_per_batch):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunks_per_batch + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id

                    current_tensor = x[bin_id].root.input_matrix[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    if param.no_indel:
                        current_tensor = np.concatenate([current_tensor[:,:,:4], current_tensor[:,:, 9:13], current_tensor[:,:,34:38], current_tensor[:,:,43:47]], axis=-1)

                    input_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = current_tensor

                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size, :param.label_size]

                # input_matrix = np.concatenate((normal_matrix, tumor_matrix), axis=1)
                # if add_contrastive:
                # label_for_normal = [[0,1] if np.argmax(item) == 1 else [1,0] for item in label]
                if discard_germline:
                    label_for_tumor = [[0,1] if np.argmax(item[:3]) == 2 else ([1,0] if np.argmax(item[:3]) == 1 else [1,0]) for item in label]
                else:
                    label_for_tumor = [[0,0,1] if np.argmax(item[:3]) == 2 else ([0,1,0] if np.argmax(item[:3]) == 1 else [1,0,0]) for item in label]
                # label_for_normal = np.array(label_for_normal, dtype=np.float32)
                label_for_tumor = np.array(label_for_tumor, dtype=np.float32)
                af_tensor = []
                if param.add_af_in_label:
                    af_list = [item[3] for item in label]
                    af_list = [[item, item, 1/float(max(item, 0.05)) / 20.0 ] for item in af_list]
                    af_tensor = torch.from_numpy(np.array(af_list)).to(device)

                # input_matrix = np.maximum(input_matrix * 2.55, 255.0)
                if not args.pileup:
                    input_tensor = torch.from_numpy(np.transpose(input_matrix, (0,3,1,2))/100.0).to(device)
                else:
                    input_tensor = torch.from_numpy(input_matrix).to(device)
                # tumor_tensor = torch.from_numpy(tumor_matrix)
                label_tensor = torch.from_numpy(label_for_tumor).to(device)
                yield input_tensor, label_tensor, af_tensor

    if args.pileup:
        from torchinfo import summary
        print (summary(model, input_size=tuple([100] + tensor_shape), device=device))
    else:
        from torchsummary import summary
        print(summary(model, input_size=(param.channel_size,param.max_depth,param.no_of_positions), device=device))
    train_dataset_loder = DataGenerator(table_dataset_list, train_shuffle_chunk_list, True)
    validate_dataset_loder = DataGenerator(validate_table_dataset_list if validation_fn else table_dataset_list, validate_shuffle_chunk_list, False)

    criterion = FocalLoss() if apply_focal_loss else nn.CrossEntropyLoss()
    criterion = criterion.to(device)
    if param.add_af_in_label:
        af_loss = AFLoss().to(device)
    # optimizer
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=param.weight_decay)
    # learning rate scheduler
    # lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min')
    lr_scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9)

    train_steps = train_data_size // batch_size
    validate_steps = validate_data_size // batch_size
    print("[INFO] Using GPU for model training: {}".format(True if device == 'cuda' else False))
    print("[INFO] The size of dataset: {}".format(train_data_size))
    print("[INFO] The training batch size: {}".format(batch_size))
    print("[INFO] The training learning_rate: {}".format(learning_rate))
    print("[INFO] The output model folder: {}".format(ochk_prefix))
    print("[INFO] Apply focal loss in training: {}".format(apply_focal_loss))
    print("[INFO] Discard germline in training: {}".format(discard_germline))
    print("[INFO] Add L2 regularization to model parameters: {}".format(add_l2_regulation_loss))
    print('[INFO] Train steps:{}'.format(train_steps))
    print('[INFO] Validate steps:{}'.format(validate_steps))

    training_loss, validation_loss = 0.0, 0.0
    training_step, validation_step = 0, 0
    echo_each_step = 200
    for epoch in range(1, max_epoch+1):
        epoch_loss = 0
        fp, tp, fn = 0,0,0
        t = tqdm(enumerate(train_dataset_loder), total=train_steps, position=0, leave=True)
        v = tqdm(enumerate(validate_dataset_loder), total=validate_steps, position=0, leave=True)
        for batch_idx, (data, label, af_list) in t:
            t.set_description('EPOCH {}'.format(epoch))
            data = data.to(device)
            label = label.to(device)

            output_logit = model(data).contiguous()
            y_truth = torch.argmax(label, axis=1)
            optimizer.zero_grad()
            loss = criterion(input=output_logit, target=label) if apply_focal_loss else criterion(output_logit, y_truth)
            # loss_2 = nn.CrossEntropyLoss().to(device)(output_logit, y_truth)
            # loss_1 = nn.CrossEntropyLoss().to(device)(output_logit, y_truth)
            # loss_2 = FocalLoss(gamma=0).to(device)(input=output_logit, target=label)
            # a = float(loss_1.detach().cpu().numpy())
            # b = float(loss_2.detach().cpu().numpy())
            # print(" ", a, b, round(a, 3) == round(b, 3))

            l2_regularization_loss = sum([l2_regularization_lambda * 0.5 * params.norm(2) ** 2 for params in model.parameters()]) if add_l2_regulation_loss else 0.0

            loss += l2_regularization_loss

            if param.add_af_in_label:
                loss += af_loss(input=output_logit, target=af_list)
            loss.backward()
            # torch.nn.utils.clip_grad_norm_(model.parameters(), param.grad_norm_clip) # apply gradient normalization
            optimizer.step()

            training_step += 1
            training_loss += loss.item()
            if add_writer:
                if training_step % echo_each_step == echo_each_step - 1:
                    writer.add_scalar('training loss', training_loss / echo_each_step, training_step)
                training_loss = 0.0
            y_truth = y_truth.cpu().numpy()
            y_pred = output_logit.argmax(dim=1).cpu().numpy()
            arg_index = 1 if discard_germline else 2
            fp += sum([True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)])
            fn += sum([True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)])
            tp += sum([True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)])

            if batch_idx + 1 == train_steps:
                break

            epoch_loss += loss
            el = epoch_loss.detach().cpu().numpy()
            t.set_postfix({'loss': el, 'tp': tp, 'fp': fp, 'fn': fn})
            t.update(1)

        # validation
        val_fp, val_tp, val_fn = 0, 0, 0
        val_epoch_loss = 0
        for batch_idx, (data, label) in v:
            v.set_description('VAL EPOCH {}'.format(epoch))
            data = data.to(device)
            label = label.to(device)
            with torch.no_grad():
                output_logit = model(data)
            y_truth = torch.argmax(label, axis=1)
            optimizer.zero_grad()
            loss = criterion(input=output_logit, target=label) if apply_focal_loss else criterion(output_logit, y_truth)
            validation_step += 1
            validation_loss += loss.item()
            if add_writer:
                if validation_step % echo_each_step == echo_each_step - 1:
                    writer.add_scalar('validation loss', validation_loss / echo_each_step, validation_step)
                validation_loss = 0.0
            y_truth = y_truth.cpu().numpy()
            y_pred = output_logit.argmax(dim=1).cpu().numpy()
            arg_index = 1 if discard_germline else 2
            val_fp += sum([True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)])
            val_fn += sum([True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)])
            val_tp += sum([True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)])

            if batch_idx + 1 == validate_steps:
                break

            val_epoch_loss += loss
            el = val_epoch_loss.detach().cpu().numpy()
            v.set_postfix(
                {'validation_loss': el, 'val_tp': val_tp, 'val_fp': val_fp, 'val_fn': val_fn})

            v.update(1)

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

        # leanrning rate decay in each epoch end
        lr_scheduler.step()
        save_path = os.path.join(ochk_prefix, "{}.pkl".format(epoch)) if ochk_prefix is not None else "{}.pkl".format(epoch)
        torch.save(model, save_path)
        # model = torch.load("{}.pkl".format(epoch))


    if add_writer:
        writer.close()
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

    parser.add_argument('--learning_rate', type=float, default=None,
                        help="Set the initial learning rate, default: %(default)s")

    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")

    parser.add_argument('--use_siam', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--add_contrastive', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="Define training samples to be excluded")
    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--use_resnet', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--add_writer', type=str2bool, default=False,
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
