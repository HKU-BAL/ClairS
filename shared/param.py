# parameters
caller_name = "clair-somatic"
version = "0.1"

from itertools import accumulate

zstd = 'gzip'

clair3_option = {
    'min_coverage': 8,
    'snp_min_af': 0.15,
    'indel_min_af': 0.2,
    'longphase_for_phasing': True
}

min_mq = ont_min_bq = 20
min_bq = 0
min_coverage = 4
split_bed_size = 10000
snv_min_af = 0.05
normal_snv_max_af = 0.05
tensor_max_depth = 168
center_padding_depth = 2
min_rescale_cov = 50
min_qual = 0.8
SAMTOOLS_VIEW_FILTER_FLAG = 2316
extend_bp = 100
alternative_base_num = min_tumor_support_read_num = 3
max_depth = tensor_max_depth + center_padding_depth
normal_tumor_ratio = 1
normal_pro = normal_tumor_ratio / (1 + normal_tumor_ratio)
max_normal_depth = int(tensor_max_depth * normal_pro)
max_tumor_depth = tensor_max_depth - max_normal_depth
ont_tensor_max_depth = 128
ont_max_depth = ont_tensor_max_depth + center_padding_depth
ont_normal_tumor_ratio = 0.7
normal_pro = ont_normal_tumor_ratio / (1 + ont_normal_tumor_ratio)
min_bq_dict = {'ont': ont_min_bq, 'ilmn': min_bq}
min_thred_qual = {'ont': 8, 'ilmn': 2}
ont_max_normal_depth = int(ont_tensor_max_depth * normal_pro)
ont_max_tumor_depth = ont_tensor_max_depth - ont_max_normal_depth

support_platform = {'ont', 'ilmn'}
matrix_depth_dict = {'ont': ont_max_depth, 'ilmn': max_depth}
normal_matrix_depth_dict = {'ont': ont_max_normal_depth, 'ilmn': max_normal_depth}
tumor_matrix_depth_dict = {'ont': ont_max_tumor_depth, 'ilmn': max_tumor_depth}
phase_normal = {'ont': False, 'ilmn': False}
phase_tumor = {'ont': True, 'ilmn': False}
qual_dict = {'ont': min_qual, 'ilmn': min_qual}
af_dict = {'ont': 0.5, 'ilmn': None}


pileup_channel = ['A', 'C', 'G', 'T', 'I', 'I1', 'D', 'D1', '*', 'a', 'c', 'g', 't', 'i', 'i1', 'd', 'd1', '#']
pileup_channel += [
    'ALMQ', 'CLMQ', 'GLMQ', 'TLMQ', 'aLMQ', 'cLMQ', 'gLMQ', 'tLMQ', 'ALBQ', 'CLBQ', 'GLBQ', 'TLBQ', 'aLBQ', 'cLBQ',
    'gLBQ', 'tLBQ']

# Full alignment input feature list
channel = [
    'reference_base', 'alternative_base', 'base_quality', 'strand_info', 'insert_base', 'phasing_info',
    "mapping_quality"]
channel_size = len(channel)
pileup_channel_size = len(pileup_channel)
tumor_channel_size = pileup_channel_size
flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
input_shape = [matrix_depth_dict['ilmn'], no_of_positions, channel_size]
ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, channel_size]


# Training hyper parameters
use_alt_base = True
label_shape = [3]
label_size = sum(label_shape)
apply_focal_loss = True
discard_germline = False
add_af_in_label = False
add_l2_regulation_loss = True
use_tf = False
smoothing = None
somatic_arg_index = 1 if discard_germline else 2
# label_size = 2
label_shape_cum = list(accumulate(label_shape))
expandReferenceRegion = 1000
NORMALIZE_NUM = 100
chunk_size = 100
trainBatchSize = 800
predictBatchSize = 250
test_chunk_size = predictBatchSize
initialLearningRate = 5e-4
l2_regularization_lambda = 1e-4
trainingDatasetPercentage = 0.8
weight_decay = 1e-6
maxEpoch = 20
OPERATION_SEED = None
RANDOM_SEED = None
NORMAL_PREFIX = 'n'
TUMOR_PREFIX = 't'
variant_type = {'ref', 'homo_somatic', 'homo_germline', 'hetero_germline'}
grad_norm_clip = 1.0

use_beta_subsampling = True
use_exp_subsampling = False

## discretized Beta cumulative distribution using numpy
beta_acc_per = [
    0.00000, 0.00288, 0.00842, 0.01639, 0.02658, 0.03880, 0.05285, 0.06856, 0.08576, 0.10427,
    0.12396, 0.14466, 0.16625, 0.18860, 0.21157, 0.23506, 0.25896, 0.28317, 0.30758, 0.33212,
    0.35669, 0.38123, 0.40566, 0.42992, 0.45394, 0.47767, 0.50105, 0.52406, 0.54663, 0.56873,
    0.59034, 0.61142, 0.63195, 0.65189, 0.67125, 0.68999, 0.70810, 0.72559, 0.74243, 0.75863,
    0.77418, 0.78908, 0.80333, 0.81695, 0.82993, 0.84228, 0.85401, 0.86513, 0.87566, 0.88560,
    0.89497, 0.90379, 0.91207, 0.91982, 0.92707, 0.93384, 0.94013, 0.94598, 0.95139, 0.95639,
    0.96099, 0.96522, 0.96910, 0.97264, 0.97586, 0.97879, 0.98143, 0.98381, 0.98595, 0.98786,
    0.98956, 0.99107, 0.99239, 0.99355, 0.99457, 0.99545, 0.99620, 0.99685, 0.99740, 0.99786,
    0.99824, 0.99855, 0.99881, 0.99902, 0.99918, 0.99931, 0.99941, 0.99949, 0.99954, 0.99958,
    0.99961, 0.99963, 0.99964, 0.99964, 0.99965, 0.99965, 0.99965, 0.99965, 0.99965, 1.00000,
]