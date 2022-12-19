import logging
import tables
import sys
import os
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

from tqdm import tqdm
from subprocess import run
from argparse import ArgumentParser, SUPPRESS
from torch.optim.lr_scheduler import StepLR
from torch.utils.tensorboard import SummaryWriter
from torch.autograd import Variable

from shared.utils import str2bool
import shared.param as param
import clair_somatic.model as model_path

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


class FocalLoss(nn.Module):
    def __init__(self, alpha=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma

    def forward(self, input, target):
        y_pred, y_true = input, target
        y_pred = torch.nn.functional.softmax(y_pred, dim=1)
        y_pred = torch.clamp(y_pred, min=1e-9, max=1 - 1e-9)
        cross_entropy = -y_true * torch.log(y_pred)
        weight = ((1 - y_pred) ** self.gamma) * y_true
        FCLoss = cross_entropy * weight

        reduce_fl = torch.mean(torch.sum(FCLoss, dim=1))
        return reduce_fl


class AFLoss(nn.Module):
    def __init__(self, alpha=None, gamma=2):
        super(AFLoss, self).__init__()
        self.gamma = gamma

    def forward(self, input, target):
        y_pred, y_true = input, target
        y_pred = torch.nn.functional.softmax(y_pred, dim=1)
        y_pred = torch.clamp(y_pred, min=1e-9, max=1 - 1e-9)

        cross_entropy = -y_true * torch.log(y_pred)
        reduce_fl = torch.mean(torch.sum(cross_entropy, dim=1))

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
    for ctg_name in ctg_name_list:
        if ctg_name + '.' in fn:
            return True
    return False


def train_model(args):
    apply_focal_loss = param.apply_focal_loss
    discard_germline = param.discard_germline
    add_l2_regulation_loss = param.add_l2_regulation_loss
    l2_regularization_lambda = param.l2_regularization_lambda
    debug_mode = args.debug_mode
    use_resnet = args.use_resnet
    platform = args.platform
    ctg_name_string = args.ctg_name
    chkpnt_fn = args.chkpnt_fn
    ochk_prefix = args.ochk_prefix
    add_writer = args.add_writer
    smoothing = args.smoothing
    phase_tumor = args.phase_tumor and args.platform == 'ont'

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
        device = 'cuda'
    apply_softmax = False if apply_focal_loss else True
    if args.pileup:
        channel_size = param.pileup_channel_size
        tumor_channel_size = param.tumor_channel_size if phase_tumor else channel_size
        pileup_tensor_shape = [param.no_of_positions, channel_size + tumor_channel_size]  # normal and tumor
        tensor_shape = pileup_tensor_shape
        model = model_path.bigru(apply_softmax=apply_softmax,
                                 num_classes=2 if discard_germline else 3,
                                 channel_size=tensor_shape[1]).to(device)

        input = torch.ones(size=[100] + tensor_shape).to(device)

    else:
        model = model_path.ResNet(platform=platform).to(device)

    if chkpnt_fn is not None:
        model = torch.load(chkpnt_fn)

    output = model(input)

    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    assert batch_size % chunk_size == 0
    chunks_per_batch = batch_size // chunk_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.max_epoch if args.max_epoch else param.maxEpoch
    bin_list = os.listdir(args.bin_fn)

    bin_list = [f for f in bin_list if
                pass_chr(f, ctg_name_list) and not exist_file_prefix(exclude_training_samples, f)]
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
            validate_chunk_num = int(
                max(1., np.floor(total_batches * (1 - training_dataset_percentage))) * chunks_per_batch)
            train_chunk_num = int(total_chunks - validate_chunk_num)
        else:
            train_chunk_num = total_chunks
        train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_chunk_num,
                                                                               chunks_per_batch,
                                                                               training_dataset_percentage)
        train_chunk_num = len(train_shuffle_chunk_list)
        validate_chunk_num = len(validate_shuffle_chunk_list)
    train_data_size = train_chunk_num * chunk_size
    validate_data_size = validate_chunk_num * chunk_size

    def DataGenerator(x, shuffle_chunk_list, train_flag=True):

        batch_num = len(shuffle_chunk_list) // chunks_per_batch
        input_matrix = np.empty([batch_size] + tensor_shape, np.float32)
        label = np.empty((batch_size, param.label_size), np.float32)
        if debug_mode:
            position_info = np.empty((batch_size, 1), "S100")
            normal_info = np.empty((batch_size, 1), "S1000")
            tumor_info = np.empty((batch_size, 1), "S1000")
        for epoch in range(max_epoch):
            random_start_position = np.random.randint(0, batch_size) if train_flag else 0
            if train_flag:
                np.random.shuffle(shuffle_chunk_list)
            for batch_idx in range(batch_num):

                for chunk_idx in range(chunks_per_batch):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunks_per_batch + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id

                    current_tensor = x[bin_id].root.input_matrix[
                                     random_start_position + chunk_id * chunk_size:random_start_position + (
                                                 chunk_id + 1) * chunk_size]
                    if not args.pileup:
                        current_tensor = current_tensor[:, :, :, :param.channel_size]

                    if param.no_indel:
                        current_tensor = np.concatenate(
                            [current_tensor[:, :, :4], current_tensor[:, :, 9:13], current_tensor[:, :, 34:38],
                             current_tensor[:, :, 43:47]], axis=-1)

                    input_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = current_tensor

                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                                                                                 random_start_position + chunk_id * chunk_size:random_start_position + (
                                                                                             chunk_id + 1) * chunk_size,
                                                                                 :param.label_size]
                    if debug_mode:
                        position_info[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position[
                                                                                             random_start_position + chunk_id * chunk_size:random_start_position + (
                                                                                                         chunk_id + 1) * chunk_size]
                        normal_info[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[
                                                                                               bin_id].root.normal_alt_info[
                                                                                           random_start_position + chunk_id * chunk_size:random_start_position + (
                                                                                                       chunk_id + 1) * chunk_size]
                        tumor_info[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.tumor_alt_info[
                                                                                          random_start_position + chunk_id * chunk_size:random_start_position + (
                                                                                                      chunk_id + 1) * chunk_size]

                positive = 1 - smoothing if smoothing is not None else 1
                negative = smoothing if smoothing is not None else 0
                if discard_germline:
                    label_for_tumor = [[negative, positive] if np.argmax(item[:3]) == 2 else \
                                           ([positive, negative] if np.argmax(item[:3]) == 1 else [positive, negative])
                                       for item in label]
                else:
                    label_for_tumor = [[negative, negative, positive] if np.argmax(item[:3]) == 2 else \
                                           ([negative, positive, negative] if np.argmax(item[:3]) == 1 else [positive,
                                                                                                             negative,
                                                                                                             negative])
                                       for item in label]

                label_for_tumor = np.array(label_for_tumor, dtype=np.float32)
                af_tensor = []
                if param.add_af_in_label:
                    af_list = [item[3] for item in label]
                    af_list = [[item, item, 1 / float(max(item, 0.05)) / 20.0] for item in af_list]
                    af_tensor = torch.from_numpy(np.array(af_list)).to(device)

                if not args.pileup:
                    input_tensor = torch.from_numpy(np.transpose(input_matrix, (0, 3, 1, 2)) / 100.0).to(device)
                else:
                    input_tensor = torch.from_numpy(input_matrix).to(device)
                label_tensor = torch.from_numpy(label_for_tumor).to(device)
                if debug_mode:
                    yield input_tensor, label_tensor, position_info, normal_info, tumor_info
                else:
                    yield input_tensor, label_tensor, af_tensor, None, None

    if args.pileup:
        try:
            from torchinfo import summary
            print(summary(model, input_size=tuple([100] + tensor_shape), device=device))
        except:
            pass
    else:
        try:
            from torchsummary import summary
            print(summary(model,
                          input_size=(param.channel_size, param.matrix_depth_dict[platform], param.no_of_positions),
                          device=device))
        except:
            pass

    train_dataset_loder = DataGenerator(table_dataset_list, train_shuffle_chunk_list, True)
    validate_dataset_loder = DataGenerator(validate_table_dataset_list if validation_fn else table_dataset_list,
                                           validate_shuffle_chunk_list, False)

    criterion = FocalLoss() if apply_focal_loss else nn.CrossEntropyLoss()
    criterion = criterion.to(device)
    if param.add_af_in_label:
        af_loss = AFLoss().to(device)

    # optimizer
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=param.weight_decay)
    # learning rate scheduler
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
    for epoch in range(1, max_epoch + 1):
        epoch_loss = 0
        fp, tp, fn = 0, 0, 0
        t = tqdm(enumerate(train_dataset_loder), total=train_steps, position=0, leave=True)
        v = tqdm(enumerate(validate_dataset_loder), total=validate_steps, position=0,
                 leave=True) if not debug_mode else enumerate(validate_dataset_loder)
        model.train()
        for batch_idx, (data, label, af_list, _, _) in t:
            t.set_description('EPOCH {}'.format(epoch))
            data = data.to(device)
            label = label.to(device)
            output_logit = model(data).contiguous()
            y_truth = torch.argmax(label, axis=1)
            optimizer.zero_grad()
            loss = criterion(input=output_logit, target=label) if apply_focal_loss else criterion(output_logit, y_truth)

            l2_regularization_loss = sum([l2_regularization_lambda * 0.5 * params.norm(2) ** 2 for params in
                                          model.parameters()]) if add_l2_regulation_loss else 0.0

            loss += l2_regularization_loss

            if param.add_af_in_label:
                loss += af_loss(input=output_logit, target=af_list)
            loss.backward()
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
        model.eval()
        for batch_idx, (data, label, position_info, normal_info, tumor_info) in v:
            if not debug_mode:
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
            output_pro = torch.softmax(output_logit, dim=1)
            y_pred = output_pro.argmax(dim=1).cpu().numpy()
            arg_index = 1 if discard_germline else 2
            if debug_mode:
                if device == 'cuda':
                    data_numpy = data.cpu().numpy()
                else:
                    data_numpy = data.numpy()
                cpu_logit = output_pro.cpu().numpy()
                for idx, (x, y) in enumerate(zip(y_truth, y_pred)):
                    if x == arg_index and y != arg_index:
                        print(idx, 'FN', x, y, cpu_logit[idx], position_info[idx], normal_info[idx], tumor_info[idx])
                        a = data_numpy[idx]
                        tmp = 1
                    if x != arg_index and y == arg_index:
                        print(idx, 'FP', x, y, cpu_logit[idx], position_info[idx], normal_info[idx], tumor_info[idx])
                        a = data_numpy[idx]
                        tmp = 1

            val_fp += sum([True if x != arg_index and y == arg_index else False for x, y in zip(y_truth, y_pred)])
            val_fn += sum([True if x == arg_index and y != arg_index else False for x, y in zip(y_truth, y_pred)])
            val_tp += sum([True if x == y and x == arg_index else False for x, y in zip(y_truth, y_pred)])

            if batch_idx + 1 == validate_steps:
                break

            val_epoch_loss += loss
            el = val_epoch_loss.detach().cpu().numpy()
            if not debug_mode:
                v.set_postfix(
                    {'validation_loss': el, 'val_tp': val_tp, 'val_fp': val_fp, 'val_fn': val_fn})

                v.update(1)

        # leanrning rate decay in each epoch end
        lr_scheduler.step()
        save_path = os.path.join(ochk_prefix, "{}.pkl".format(epoch)) if ochk_prefix is not None else "{}.pkl".format(
            epoch)
        print(save_path)
        torch.save(model, save_path)

    if add_writer:
        writer.close()
    for table_dataset in table_dataset_list:
        table_dataset.close()

    for table_dataset in validate_table_dataset_list:
        table_dataset.close()

def main():
    parser = ArgumentParser(description="Train a somatic model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input, default: %(default)s")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--bin_fn', type=str, default="", required=True,
                        help="Binary tensor input, support multiple bin readers using pytables")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a model to resume training or for fine-tuning")

    parser.add_argument('--ochk_prefix', type=str, default=None,
                        help="Prefix for model output after each epoch")

    # options for advanced users
    parser.add_argument('--max_epoch', type=int, default=None,
                        help="Maximum number of training epochs")

    parser.add_argument('--smoothing', type=int, default=param.smoothing,
                        help="Label smoothing Epsilon in training")

    parser.add_argument('--learning_rate', type=float, default=None,
                        help="Set the initial learning rate, default: %(default)s")

    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")

    # mutually-incompatible validation options
    vgrp = parser.add_mutually_exclusive_group()
    vgrp.add_argument('--random_validation', action='store_true',
                      help="Use random sample of dataset for validation, default: %(default)s")

    vgrp.add_argument('--validation_fn', type=str, default=None,
                      help="Binary tensor input for use in validation: %(default)s")

    # Internal process control
    ## use siamese network in training
    parser.add_argument('--use_siam', action='store_true',
                        help=SUPPRESS)

    ## use contrastive loss in training
    parser.add_argument('--add_contrastive', action='store_true',
                        help=SUPPRESS)

    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    ## use resnet for model training
    parser.add_argument('--use_resnet', type=str2bool, default=False,
                        help=SUPPRESS)

    ## add logging writer using torchvision
    parser.add_argument('--add_writer', type=str2bool, default=False,
                        help=SUPPRESS)

    ## Debug mode
    parser.add_argument('--debug_mode', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--phase_tumor', type=str2bool, default=False,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    train_model(args)


if __name__ == "__main__":
    main()
