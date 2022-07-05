import sys
import os
import math
import tables
# import tensorflow as tf
import numpy as np
import logging
import torch
import shlex
from time import time
from argparse import ArgumentParser, SUPPRESS
from threading import Thread
from subprocess import PIPE, run
from clair_somatic.call_variants import output_vcf_from_probability, OutputConfig
from math import log, e
from collections import namedtuple

from clair_somatic.task.gt21 import (
    GT21_Type, gt21_enum_from_label,
    HOMO_SNP_GT21, HOMO_SNP_LABELS,
    HETERO_SNP_GT21, HETERO_SNP_LABELS, GT21_LABELS, partial_label_from, mix_two_partial_labels
)
import clair_somatic.utils as utils
from clair_somatic.task.genotype import Genotype, genotype_string_from, genotype_enum_from, genotype_enum_for_task
from shared.utils import IUPAC_base_to_ACGT_base_dict as BASE2ACGT, BASIC_BASES, str2bool, file_path_from, log_error, log_warning, subprocess_popen, TensorStdout
from clair_somatic.task.variant_length import VariantLength
import shared.param as param
# os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

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

def print_output_message(
            output_file,
            chromosome,
            position,
            reference_base,
            normal_alt_info,
            tumor_alt_info,
            gt21_probabilities,
            genotype_probabilities,
            variant_length_probabilities_1,
            variant_length_probabilities_2,
            extra_infomation_string=""
    ):
        global call_fn
        if call_fn is not None:
            output_vcf_from_probability(
            chromosome,
            position,
            reference_base,
            normal_alt_info,
            tumor_alt_info,
            gt21_probabilities,
            output_config=output_config,
            vcf_writer=output_file
            )
        else:
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                chromosome,
                position,
                reference_base,
                normal_alt_info,
                tumor_alt_info,
                ' '.join(["{:0.8f}".format(x) for x in gt21_probabilities]),
                # ["{:0.8f}".format(x) for x in genotype_probabilities],
                # ["{:0.8f}".format(x) for x in variant_length_probabilities_1],
                # ["{:0.8f}".format(x) for x in variant_length_probabilities_2],
                extra_infomation_string
            ))

def tensor_generator_from(tensor_file_path, batch_size, pileup=False, platform='ont'):
    float_type = 'float32'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    if pileup:
        channel_size = param.pileup_channel_size
        tensor_shape = [param.no_of_positions, channel_size * 2]
    else:
        tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        contig, coord, seq, normal_tensor, normal_alt_info, tumor_tensor, tumor_alt_info, variant_type = row.split("\t")
        # if pileup:
        #     tensor = np.array(tensor.split(), dtype=np.dtype(float_type))
        #     depth = int(alt_info.split('-', maxsplit=1)[0])
        #     max_depth = param.max_depth_dict[platform]
        #     # for extreme high coverage data, make sure we could have a truncated coverage
        #     if depth > 0 and depth > max_depth * 1.5:
        #         scale_factor = depth / max_depth
        #         tensor = tensor / scale_factor
        # else:
            # need add padding if depth is lower than maximum depth.
        normal_matrix = [float(item) for item in normal_tensor.split()]
        tumor_matrix = [float(item) for item in tumor_tensor.split()]

        if pileup:
            apply_normalize = False
            channel_size = param.pileup_channel_size
            tensor = []
            for idx in range(param.no_of_positions):
                if apply_normalize:
                    normal_coverage = float(normal_alt_info.split('-')[0])
                    tumor_coverage = float(tumor_alt_info.split('-')[0])
                    tensor += [float(item) / normal_coverage for item in
                                     normal_matrix[idx * channel_size: (idx + 1) * channel_size]]
                    tensor += [float(item) / tumor_coverage for item in
                                     tumor_matrix[idx * channel_size: (idx + 1) * channel_size]]
                else:
                    tensor += normal_matrix[idx * channel_size: (idx + 1) * channel_size]
                    tensor += tumor_matrix[idx * channel_size: (idx + 1) * channel_size]
        else:
            normal_depth = len(normal_matrix) // tensor_shape[1] // tensor_shape[2]
            tumor_depth = len(tumor_matrix) // tensor_shape[1] // tensor_shape[2]
            tensor_depth = normal_depth + tumor_depth
            center_padding_depth = param.center_padding_depth
            padding_depth = tensor_shape[0] - tensor_depth - center_padding_depth
            prefix_padding_depth = int(padding_depth / 2)
            suffix_padding_depth = padding_depth - int(padding_depth / 2)
            prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
            center_zero_padding = [0] * center_padding_depth * tensor_shape[1] * tensor_shape[2]
            suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
            tensor = prefix_zero_padding + normal_matrix + center_zero_padding + tumor_matrix + suffix_zero_padding
        tensor = np.array(tensor, dtype=np.dtype(float_type))

        # tensor_depth = len(tensor) // tensor_shape[1] // tensor_shape[2]
        # padding_depth = tensor_shape[0] - tensor_depth
        # prefix_padding_depth = int(padding_depth / 2)
        # suffix_padding_depth = padding_depth - int(padding_depth / 2)
        # prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
        # suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
        # tensor = prefix_zero_padding + tensor + suffix_zero_padding
        # tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = contig + ":" + coord + ":" + seq
        return tensor, pos, seq, normal_alt_info, tumor_alt_info, variant_type

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        normal_alt_info_list = []
        tumor_alt_info_list = []
        variant_type_list = []
        for tensor, pos, seq, normal_alt_info, tumor_alt_info, variant_type in batch:
            if seq[param.flankingBaseNum] not in "ACGT":
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            normal_alt_info_list.append(normal_alt_info)
            tumor_alt_info_list.append(tumor_alt_info)
            variant_type_list.append(variant_type)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= 0:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], normal_alt_info_list[:current_batch_size],tumor_alt_info_list[:current_batch_size], variant_type_list[:current_batch_size]
        # for p, (pos, n_info, t_info) in enumerate(
        #         zip(positions[:current_batch_size], normal_alt_info_list[:current_batch_size],
        #             tumor_alt_info_list[:current_batch_size])):
        #     print(p, pos, n_info, t_info)
    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


# def Run(args):
#     os.environ["OMP_NUM_THREADS"] = "1"
#     os.environ["OPENBLAS_NUM_THREADS"] = "1"
#     os.environ["MKL_NUM_THREADS"] = "1"
#     os.environ["NUMEXPR_NUM_THREADS"] = "1"
#
#     # tf.config.threading.set_intra_op_parallelism_threads(1)
#     # tf.config.threading.set_inter_op_parallelism_threads(1)
#
#     global test_pos
#     test_pos = None
#
#     predict(args=args)

def batch_output_for_ensemble(X, batch_chr_pos_seq, alt_info_list, batch_Y, output_config, output_utilities):
    batch_size = len(batch_chr_pos_seq)
    batch_gt21_probabilities, batch_genotype_probabilities, = batch_Y

    if len(batch_gt21_probabilities) != batch_size:
        sys.exit(
            "Inconsistent shape between input tensor and output predictions %d/%d" %
            (batch_size, len(batch_gt21_probabilities))
        )

    tensor_position_center = param.flankingBaseNum

    for (
            x,
            chr_pos_seq,
            gt21_probabilities,
            genotype_probabilities,
            alt_info
    ) in zip(
        X,
        batch_chr_pos_seq,
        batch_gt21_probabilities,
        batch_genotype_probabilities,
        alt_info_list
    ):
        if output_config.tensor_fn != 'PIPE':
            chromosome, position, reference_sequence = chr_pos_seq.decode().rstrip().split(":")
        else:
            chromosome, position, reference_sequence = chr_pos_seq

        position = int(position)

        if reference_sequence[tensor_position_center] not in BASIC_BASES:
            continue

        output_utilities.output(
            "\t".join(
                [
                    chromosome,
                    str(position),
                    reference_sequence,
                    alt_info.decode(),
                    ' '.join(["{:0.6f}".format(p) for p in list(gt21_probabilities)]),
                    ' '.join(["{:0.6f}".format(p) for p in list(genotype_probabilities)])]
            )
        )


def batch_output(output_file, batch_chr_pos_seq, normal_alt_info_list, tumor_alt_info_list, batch_Y):
    batch_size = len(batch_chr_pos_seq)

    batch_gt21_probabilities = batch_Y[:,:param.label_shape_cum[0]]
    # batch_genotype_probabilities = batch_Y[:,param.label_shape_cum[0]:param.label_shape_cum[1]]
    batch_genotype_probabilities = [0] * batch_size
    if len(batch_gt21_probabilities) != batch_size:
        sys.exit(
            "Inconsistent shape between input tensor and output predictions %d/%d" %
            (batch_size, len(batch_gt21_probabilities))
        )
    batch_variant_length_probabilities_1, batch_variant_length_probabilities_2 = [0] * batch_size, [0] * batch_size

    # if output_config.add_indel_length:
    #     batch_variant_length_probabilities_1, batch_variant_length_probabilities_2 = batch_Y[:,param.label_shape_cum[1]:param.label_shape_cum[2]], batch_Y[:,param.label_shape_cum[2]:param.label_shape_cum[3]]
    for (
            chr_pos_seq,
            normal_alt_info,
            tumor_alt_info,
            gt21_probabilities,
            genotype_probabilities,
            variant_length_probabilities_1,
            variant_length_probabilities_2
    ) in zip(
        batch_chr_pos_seq,
        normal_alt_info_list,
        tumor_alt_info_list,
        batch_gt21_probabilities,
        batch_genotype_probabilities,
        batch_variant_length_probabilities_1,
        batch_variant_length_probabilities_2
    ):

        output_with(
            output_file,
            chr_pos_seq,
            normal_alt_info,
            tumor_alt_info,
            gt21_probabilities,
            genotype_probabilities,
            variant_length_probabilities_1,
            variant_length_probabilities_2,
        )


def output_with(
        output_file,
        chr_pos_seq,
        normal_alt_info,
        tumor_alt_info,
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
):
    if type(chr_pos_seq) == np.ndarray:
        chr_pos_seq = chr_pos_seq[0].decode()
        normal_alt_info = normal_alt_info[0].decode()
        tumor_alt_info = tumor_alt_info[0].decode()
    elif type(chr_pos_seq) == np.bytes_ or type(chr_pos_seq) == bytes:
        chr_pos_seq = chr_pos_seq.decode()
        normal_alt_info = normal_alt_info.decode()
        tumor_alt_info = tumor_alt_info.decode()


    chromosome, position, reference_sequence = chr_pos_seq.rstrip().split(':')[:3]
    reference_base = reference_sequence[param.flankingBaseNum].upper()
    print_output_message(
        output_file,
        chromosome,
        position,
        reference_base,
        normal_alt_info,
        tumor_alt_info,
        gt21_probabilities,
        genotype_probabilities,
        variant_length_probabilities_1,
        variant_length_probabilities_2,
        ""
    )

def DataGenerator(dataset, num_epoch, batch_size, chunk_start_pos, chunk_end_pos):
    for idx in range(num_epoch):
        start_pos = chunk_start_pos + idx * batch_size
        end_pos = min(chunk_start_pos + (idx + 1) * batch_size, chunk_end_pos)
        input_matrix = dataset.input_matrix[start_pos:end_pos]
        position = dataset.position[start_pos:end_pos]  # .flatten()
        normal_alt_info_list = dataset.normal_alt_info[start_pos:end_pos]  # .flatten()
        tumor_alt_info_list = dataset.tumor_alt_info[start_pos:end_pos]  # .flatten()
        yield input_matrix, position, normal_alt_info_list, tumor_alt_info_list


def predict(args):
    global output_config
    global call_fn
    output_config = OutputConfig(
        is_show_reference=args.show_ref,
        is_show_germline=args.show_germline,
        is_haploid_precise_mode_enabled=args.haploid_precise,
        is_haploid_sensitive_mode_enabled=args.haploid_sensitive,
        is_output_for_ensemble=args.output_for_ensemble,
        quality_score_for_pass=args.qual,
        tensor_fn=args.tensor_fn,
        input_probabilities=args.input_probabilities,
        add_indel_length=args.add_indel_length,
        gvcf=args.gvcf,
        pileup=args.pileup
    )

    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    predict_fn = args.predict_fn
    use_gpu = args.use_gpu
    logging.info("[INFO] Make prediction")
    variant_call_start_time = time()
    call_fn = args.call_fn
    chkpnt_fn = args.chkpnt_fn
    tensor_fn = args.tensor_fn
    platform = args.platform
    torch.set_num_threads(1)
    if use_gpu and not torch.cuda.is_available():
        print("[WARNING] --use_gpu is enabled, but cuda is not found")
        use_gpu = False
    if use_gpu:
        device = 'cuda'
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        device = 'cpu'


    if call_fn is not None:
        from shared.vcf import VcfWriter
        call_dir = os.path.dirname(call_fn)
        if not os.path.exists(call_dir):
            output = run("mkdir -p {}".format(call_dir), shell=True)
        vcf_writer = VcfWriter(vcf_fn=args.call_fn,
                               ref_fn=args.ref_fn,
                               ctg_name=args.ctg_name,
                               show_ref_calls=args.show_ref,
                               sample_name=args.sample_name,
                               )
        output_file = vcf_writer
    elif predict_fn != "PIPE":
        predict_dir = os.path.dirname(predict_fn)
        if not os.path.exists(predict_dir):
            output = run("mkdir -p {}".format(predict_dir), shell=True)
        predict_fn_fpo = open(predict_fn, "wb")
        predict_fn_fp = subprocess_popen(shlex.split("{} -c".format(param.zstd)), stdin=PIPE, stdout=predict_fn_fpo)
        output_file = predict_fn_fp.stdin
    else:
        predict_fn_fp = TensorStdout(sys.stdout)
        output_file = predict_fn_fp.stdin

    global test_pos
    test_pos = None

    # batch_output_method = batch_output_for_ensemble if output_config.is_output_for_ensemble else batch_output
    batch_output_method = batch_output
    # m.load_weights(args.chkpnt_fn)
    model = torch.load(chkpnt_fn, map_location=torch.device(device))
    total = 0
    softmax = torch.nn.Softmax(dim=1)
    if not args.is_from_tables:
        is_finish_loaded_all_mini_batches = False
        mini_batches_loaded = []
        mini_batches_to_output = []

        def load_mini_batch():
            try:
                mini_batches_loaded.append(next(tensor_generator))
            except StopIteration:
                return

        tensor_generator = tensor_generator_from(tensor_fn, param.predictBatchSize, pileup=args.pileup,
                                                       platform=platform)

        while True:
            thread_pool = []
            if len(mini_batches_to_output) > 0:
                mini_batch = mini_batches_to_output.pop(0)
                input_tensor, position, normal_alt_info_list, tumor_alt_info_list, variant_type_list = mini_batch

                if param.use_tf:
                    prediction = model.predict_on_batch(input_tensor)[0]
                else:
                    if args.pileup:
                        input_matrix = torch.from_numpy(input_tensor).to(device)
                    else:
                        input_matrix = torch.from_numpy(np.transpose(input_tensor, (0, 3, 1, 2)) / 100.0).float().to(device)
                        # print(np.unique(input_tensor[:,:,:,5]))
                    with torch.no_grad():
                        prediction = model(input_matrix)
                    prediction = softmax(prediction)
                    prediction = prediction.cpu().numpy()

                total += len(input_tensor)
                thread_pool.append(Thread(
                    target=batch_output,
                    args=(output_file, position, normal_alt_info_list, tumor_alt_info_list, prediction)
                ))

            if not is_finish_loaded_all_mini_batches:
                thread_pool.append(Thread(target=load_mini_batch))

            for t in thread_pool:
                t.start()
            for t in thread_pool:
                t.join()

            is_finish_loaded_all_mini_batches = len(mini_batches_loaded) == 0
            while len(mini_batches_loaded) > 0:
                mini_batch = mini_batches_loaded.pop(0)
                mini_batches_to_output.append(mini_batch)

            is_nothing_to_predict_and_output = (
                    len(thread_pool) <= 0 and len(mini_batches_to_output) <= 0
            )
            if is_finish_loaded_all_mini_batches and is_nothing_to_predict_and_output:
                break
    else:
        if not os.path.exists(args.tensor_fn):
            logging.info("skip {}, not existing chunk_id".format(args.tensor_fn))
            return
        dataset = tables.open_file(tensor_fn, 'r').root
        batch_size = param.predictBatchSize
        dataset_size = len(dataset.label)
        chunk_start_pos, chunk_end_pos = 0, dataset_size
        tensor_shape = param.ont_input_shape if args.platform == 'ont' else param.input_shape
        # process by chunk windows
        if chunk_id is not None and chunk_num is not None:
            chunk_dataset_size = dataset_size // chunk_num if dataset_size % chunk_num == 0 else dataset_size // chunk_num + 1
            chunk_start_pos = chunk_id * chunk_dataset_size
            dataset_size = min(chunk_dataset_size, dataset_size - chunk_start_pos)
            chunk_end_pos = min(chunk_start_pos + dataset_size, chunk_end_pos)
        num_epoch = dataset_size // batch_size if dataset_size % batch_size == 0 else dataset_size // batch_size + 1
        data_generator = DataGenerator(dataset, num_epoch, batch_size, chunk_start_pos, chunk_end_pos)

        dataset_iter = iter(data_generator)
        for idx in range(num_epoch):
            input_tensor, position, normal_alt_info_list, tumor_alt_info_list = next(dataset_iter)
            input_matrix = torch.from_numpy(np.transpose(input_tensor, (0, 3, 1, 2)) / 100.0).float().to(device)
            with torch.no_grad():
                prediction = model(input_matrix)
            prediction = softmax(prediction)
            prediction = prediction.cpu().numpy()
            start_pos = idx * batch_size
            end_pos = min((idx + 1) * batch_size, dataset_size)
            batch_output(output_file, position, normal_alt_info_list, tumor_alt_info_list, prediction)
            total += len(input_tensor)
    if predict_fn != "PIPE":
        predict_fn_fp.stdin.close()
        predict_fn_fp.wait()
        predict_fn_fpo.close()

    logging.info("Total process positions: {}".format(total))
    logging.info("Total time elapsed: %.2f s" % (time() - variant_call_start_time))


def main():
    parser = ArgumentParser(description="Call variants using a trained model and tensors of candidate variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a trained model for variant calling, required")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required if --gvcf is enabled")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--ctg_start', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--ctg_end', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctg_name if not set")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=0,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    # options for advanced users
    parser.add_argument('--temp_file_dir', type=str, default='./',
                        help="EXPERIMENTAL: The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--haploid_precise', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant")

    parser.add_argument('--haploid_sensitive', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant")

    # options for debug purpose
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building your own pipeline")

    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    parser.add_argument('--input_probabilities', action='store_true',
                        help="DEBUG: Use network probability outputs as input and generate variants from them")

    parser.add_argument('--output_probabilities', action='store_true',
                        help="DEBUG: Output the network probabilities of gt21, genotype, indel_length_1 and indel_length_2")

    # options for internal process control
    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Include indel length in training and calling, false for pileup and true for raw alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Generating outputs for ensemble model calling
    parser.add_argument('--output_for_ensemble', action='store_true',
                        help=SUPPRESS)
    
    ## Use bin file from pytables to speed up calling.
    parser.add_argument('--is_from_tables', type=str2bool, default=False,
                        help=SUPPRESS)

    ## Output reference calls
    parser.add_argument('--show_ref', action='store_true',
                        help=SUPPRESS)
    parser.add_argument('--show_germline', action='store_true',
                        help=SUPPRESS)
    args = parser.parse_args()

    # if len(sys.argv[1:]) == 0:
    #     parser.print_help()
    #     sys.exit(1)

    predict(args)


if __name__ == "__main__":
    main()

# /mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3 ${CS} predict > tmp_alt
# /autofs/bal33/zxzheng/env/miniconda2/envs/torch_3090/bin/python3 /autofs/bal36/zxzheng/somatic/Clair-somatic/clair-somatic.py predict --tensor_fn /autofs/bal36/zxzheng/somatic/ilmn/ilmn/chr20_add_germline_min_coverage_3_mq_20/build/tensor_can/chr20.13_2_3 --call_fn /autofs/bal36/zxzheng/somatic/ilmn/ilmn/chr20_add_germline_min_coverage_3_mq_20/predict/tmp.vcf --chkpnt_fn /autofs/bal36/zxzheng/somatic/ilmn/ilmn/training_add_germline_add_mq20/train/test_no_germline/12.pkl --use_gpu 1 --platform ilmn