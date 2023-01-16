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

import sys
import os
import tables
import numpy as np
import logging
import torch
import shlex

from time import time
from argparse import ArgumentParser, SUPPRESS
from threading import Thread
from sys import stderr
from subprocess import PIPE, run, Popen

from clairs.call_variants import output_vcf_from_probability, OutputConfig
from shared.utils import IUPAC_base_to_ACGT_base_dict as BASE2ACGT, BASIC_BASES, str2bool, file_path_from, log_error, \
    log_warning, subprocess_popen, TensorStdout
import shared.param as param


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
        probabilities,
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
            probabilities,
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
            ' '.join(["{:0.8f}".format(x) for x in probabilities]),
            extra_infomation_string
        ))


def tensor_generator_from(tensor_file_path, batch_size, pileup=False, min_rescale_cov=None, phase_tumor=False,
                          platform='ont'):
    float_type = 'float32'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    if pileup:
        channel_size = param.pileup_channel_size
        tumor_channel_size = param.tumor_channel_size if phase_tumor else channel_size
        tensor_shape = [param.no_of_positions, channel_size + tumor_channel_size]
    else:
        tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        contig, coord, seq, normal_tensor, normal_alt_info, tumor_tensor, tumor_alt_info, variant_type = row.split("\t")
        normal_matrix = [float(item) for item in normal_tensor.split()]
        tumor_matrix = [float(item) for item in tumor_tensor.split()]

        if pileup:
            apply_normalize = False
            if min_rescale_cov is not None:
                normal_coverage = float(normal_alt_info.split('-')[0])
                tumor_coverage = float(tumor_alt_info.split('-')[0])
                normal_rescale = float(min_rescale_cov) / normal_coverage if normal_coverage > min_rescale_cov else None
                tumor_rescale = float(min_rescale_cov) / tumor_coverage if tumor_coverage > min_rescale_cov else None

            channel_size = param.pileup_channel_size
            tumor_channel_size = param.tumor_channel_size
            tensor = []
            for idx in range(param.no_of_positions):
                if apply_normalize:

                    tensor += [float(item) / normal_coverage for item in
                               normal_matrix[idx * channel_size: (idx + 1) * channel_size]]
                    tensor += [float(item) / tumor_coverage for item in
                               tumor_matrix[idx * channel_size: (idx + 1) * channel_size]]
                else:
                    if normal_rescale is not None:
                        tensor += [item * normal_rescale for item in
                                   normal_matrix[idx * channel_size: (idx + 1) * channel_size]]
                    else:
                        tensor += normal_matrix[idx * channel_size: (idx + 1) * channel_size]

                    if tumor_rescale is not None:
                        tensor += [item * tumor_rescale for item in
                                   tumor_matrix[idx * tumor_channel_size: (idx + 1) * tumor_channel_size]]
                    else:
                        tensor += tumor_matrix[idx * tumor_channel_size: (idx + 1) * tumor_channel_size]

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
        yield X[:current_batch_size], positions[:current_batch_size], normal_alt_info_list[
                                                                      :current_batch_size], tumor_alt_info_list[
                                                                                            :current_batch_size], variant_type_list[
                                                                                                                  :current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


def get_bins(tensor_file_path, batch_size=10000, pileup=False, platform='ont'):

    def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
        return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True,
                     shell=True)

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format('gzip', tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    if pileup:
        channel_size = param.pileup_channel_size
        tensor_shape = [param.no_of_positions, channel_size * 2]
    else:
        tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    prod_tensor_shape = np.prod(tensor_shape)

    tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
    positions = []
    normal_alt_info_list = []
    tumor_alt_info_list = []
    variant_type_list = []

    for row in fo.readline():

        contig, coord, seq, normal_tensor, normal_alt_info, tumor_tensor, tumor_alt_info, variant_type = row.split("\t")
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

        pos = contig + ":" + coord + ":" + seq

        # for tensor, pos, seq, normal_alt_info, tumor_alt_info, variant_type in batch:
        if seq[param.flankingBaseNum] not in "ACGT":
            continue
        tensors[len(positions)] = tensor
        positions.append(pos)
        normal_alt_info_list.append(normal_alt_info)
        tumor_alt_info_list.append(tumor_alt_info)
        variant_type_list.append(variant_type)

    X = np.reshape(tensors, ([batch_size] + tensor_shape))

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


def batch_output(output_file, batch_chr_pos_seq, normal_alt_info_list, tumor_alt_info_list, batch_Y):
    batch_size = len(batch_chr_pos_seq)

    batch_probabilities = batch_Y[:, :param.label_shape_cum[0]]
    if len(batch_probabilities) != batch_size:
        sys.exit(
            "Inconsistent shape between input tensor and output predictions %d/%d" %
            (batch_size, len(batch_probabilities))
        )

    for (
            chr_pos_seq,
            normal_alt_info,
            tumor_alt_info,
            probabilities,
    ) in zip(
        batch_chr_pos_seq,
        normal_alt_info_list,
        tumor_alt_info_list,
        batch_probabilities,
    ):
        output_with(
            output_file,
            chr_pos_seq,
            normal_alt_info,
            tumor_alt_info,
            probabilities,
        )


def output_with(
        output_file,
        chr_pos_seq,
        normal_alt_info,
        tumor_alt_info,
        probabilities,
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
        probabilities,
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
        is_output_for_ensemble=args.output_for_ensemble,
        quality_score_for_pass=args.qual,
        tensor_fn=args.tensor_fn,
        input_probabilities=args.input_probabilities,
        pileup=args.pileup
    )

    param.flankingBaseNum = param.flankingBaseNum if args.flanking is None else args.flanking
    param.no_of_positions = param.flankingBaseNum * 2 + 1
    param.min_rescale_cov = param.min_rescale_cov if args.min_rescale_cov is None else args.min_rescale_cov
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    predict_fn = args.predict_fn
    use_gpu = args.use_gpu
    variant_call_start_time = time()
    call_fn = args.call_fn
    chkpnt_fn = args.chkpnt_fn
    tensor_fn = args.tensor_fn
    platform = args.platform
    torch.set_num_threads(1)
    torch.manual_seed(0)
    np.random.seed(0)
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

    if param.use_tf:
        import clairs.model_tf as model_path
        model = model_path.Clair3_P()
        model.load_weights(args.chkpnt_fn)

    else:
        model = torch.load(chkpnt_fn, map_location=torch.device(device))

        model.eval()

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

        tensor_generator = tensor_generator_from(tensor_file_path=tensor_fn,
                                                 batch_size=param.predictBatchSize,
                                                 pileup=args.pileup,
                                                 min_rescale_cov=param.min_rescale_cov,
                                                 phase_tumor=args.phase_tumor,
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
                        input_matrix = torch.from_numpy(np.transpose(input_tensor, (0, 3, 1, 2)) / 100.0).float().to(
                            device)
                        if input_matrix.shape[1] != param.channel_size:
                            input_matrix = input_matrix[:, :param.channel_size, :, :]
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
            batch_output(output_file, position, normal_alt_info_list, tumor_alt_info_list, prediction)
            total += len(input_tensor)

    run_time = "%.1fs" % (time() - variant_call_start_time)
    logging.info("[INFO] {} total processed positions: {}, time elapsed: {}".format(args.ctg_name, total, run_time))

    if call_fn is not None:
        output_file.close()
        if os.path.exists(call_fn):
            vcf_file = open(call_fn, 'r').readlines()
            if not len(vcf_file):
                os.remove(call_fn)
            for row in vcf_file:
                if row[0] != '#':
                    return
            logging.info("[INFO] No variant output for {}, remove empty VCF".format(call_fn.split('/')[-1]))
            os.remove(call_fn)
    elif predict_fn != "PIPE":
        predict_fn_fp.stdin.close()
        predict_fn_fp.wait()
        predict_fn_fpo.close()


def main():
    parser = ArgumentParser(description="Candidate variants probability prediction using tensors and a trained model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a trained model for calling, required")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument('--show_germline', action='store_true',
                        help="Show germline calls in VCF file")

    # options for advanced users
    parser.add_argument('--min_rescale_cov', type=int, default=param.min_rescale_cov,
                        help="EXPERIMENTAL: Minimum coverage after rescalling from excessively high coverage data")

    # options for debug purpose
    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    parser.add_argument('--input_probabilities', action='store_true',
                        help="DEBUG: Use network probability outputs as input and generate variants from them")

    parser.add_argument('--output_probabilities', action='store_true',
                        help="DEBUG: Output the network probabilities of gt21, genotype, indel_length_1 and indel_length_2")

    # options for internal process control
    ## Use GPU for calling
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help=SUPPRESS)

    ## If set, variants with >=QUAL will be marked 'PASS', or 'LowQual'
    parser.add_argument('--qual', type=int, default=0,
                        help=SUPPRESS)

    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
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

    parser.add_argument('--phase_tumor', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--flanking', type=int, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    predict(args)


if __name__ == "__main__":
    main()
