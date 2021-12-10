# Clair3 full alignment parameters
REPO_NAME = "Clair3"
from itertools import accumulate

zstd='/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/zstd'
default_optimizer = "Radam"
default_loss_function = "FocalLoss"

tensor_max_depth = 128
center_padding_depth = 2
max_depth = tensor_max_depth + center_padding_depth
normal_tumor_ratio = 0.33
min_af = 0.1
min_af_dict = {'ont':0.15, 'hifi':min_af, 'ilmn':min_af }
max_normal_depth = int(tensor_max_depth * 0.33)
max_tumor_depth = tensor_max_depth - max_normal_depth
support_platform = {'ont', 'hifi','ilmn'}
matrix_depth_dict = {'ont': max_depth, 'hifi': max_depth, 'ilmn': max_depth}
normal_matrix_depth_dict = {'ont': max_normal_depth, 'hifi': max_normal_depth, 'ilmn': max_normal_depth}
tumor_matrix_depth_dict = {'ont': max_tumor_depth, 'hifi': max_tumor_depth, 'ilmn': max_tumor_depth}
min_tumor_support_read_num = 2

# Full alignment input feature list
channel = (
'reference_base', 'alternative_base', 'base_quality', 'strand_info', 'variant_type', 'insert_base',
'phasing_info')  # phasing info if add_phasing, 'mapping_quality',
channel_size = len(channel)
flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
input_shape = [matrix_depth_dict['hifi'], no_of_positions, channel_size]
ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, channel_size]
# label_shape = [21, 3, no_of_positions, no_of_positions]
label_shape = [3]
label_size = sum(label_shape)
apply_focal_loss = True
discard_germline = True
add_l2_regulation_loss = True
somatic_arg_index = 1
# label_size = 2
label_shape_cum = list(accumulate(label_shape))
expandReferenceRegion = 1000
SAMTOOLS_VIEW_FILTER_FLAG = 2316
NORMALIZE_NUM = 100

# Realignment parameters
partition_size = 500000
realign_chunk_size = 5000
phasing_window_size = 100000
illumina_phasing_window_size = 10000
max_phasing_depth = 15
min_phasing_read_coverage = 2
split_region_size = 1000
extend_bp = 100

min_mq = 10
min_bq = 0
min_coverage = 4
split_bed_size = 10000
# Training hyperparameters
chunk_size = 50
trainBatchSize = 400
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
variant_type={'ref', 'homo_somatic', 'homo_germline', 'hete_germline'}
grad_norm_clip = 1.0