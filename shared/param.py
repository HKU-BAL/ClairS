# Clair3 full alignment parameters
REPO_NAME = "Clair3"
from itertools import accumulate
from collections import namedtuple
zstd='gzip'
default_optimizer = "Radam"
default_loss_function = "FocalLoss"

tensor_max_depth = 168
center_padding_depth = 2
max_depth = tensor_max_depth + center_padding_depth
normal_tumor_ratio = 1
normal_pro = normal_tumor_ratio/(1+normal_tumor_ratio)
min_af = 0.1
snv_min_af = 0.05
min_af_dict = {'ont':0.15, 'hifi':min_af, 'ilmn':min_af}
max_normal_depth = int(tensor_max_depth * normal_pro)
max_tumor_depth = tensor_max_depth - max_normal_depth



ont_tensor_max_depth = 128
center_padding_depth = 2
ont_max_depth = ont_tensor_max_depth + center_padding_depth
ont_normal_tumor_ratio = 0.7
normal_pro = ont_normal_tumor_ratio/(1+ont_normal_tumor_ratio)
min_af = 0.1
snv_min_af = 0.05
min_af_dict = {'ont':0.15, 'hifi':min_af, 'ilmn':min_af}
ont_max_normal_depth = int(ont_tensor_max_depth * normal_pro)
ont_max_tumor_depth = ont_tensor_max_depth - ont_max_normal_depth

support_platform = {'ont', 'hifi','ilmn'}
matrix_depth_dict = {'ont': ont_max_depth, 'hifi': max_depth, 'ilmn': max_depth}
normal_matrix_depth_dict = {'ont': ont_max_normal_depth, 'hifi': max_normal_depth, 'ilmn': max_normal_depth}
tumor_matrix_depth_dict = {'ont': ont_max_tumor_depth, 'hifi': max_tumor_depth, 'ilmn': max_tumor_depth}
min_tumor_support_read_num = 3
alternative_base_num = 3

# Full alignment input feature list
add_CP_channel = True
channel = [
'reference_base', 'alternative_base', 'base_quality', 'strand_info', 'insert_base',
'phasing_info', "mapping_quality"]  # phasing info if add_phasing, 'mapping_quality',
channel = channel + ['candidate_proportion'] if add_CP_channel else channel
channel_size = len(channel)

#                  0    1    2    3    4    5    6    7     8    9    10   11  12   13    14  15   16    17  18      19      20      21
pileup_channel = ['A', 'C', 'G', 'T', 'I', 'I1', 'D', 'D1', '*', 'a', 'c', 'g','t', 'i', 'i1','d', 'd1','#']
# channel = ['A', 'C', 'G', 'T', 'I', 'D', '*', 'a', 'c', 'g','t', 'i','d','#']
pileup_channel += [
 'ALMQ', 'CLMQ', 'GLMQ', 'TLMQ', 'aLMQ', 'cLMQ', 'gLMQ', 'tLMQ', 'ALBQ', 'CLBQ', 'GLBQ', 'TLBQ', 'aLBQ', 'cLBQ', 'gLBQ', 'tLBQ']
no_indel = False #True
no_mq = False # True
no_bq = False #True
pileup_channel_size = len(pileup_channel)
pileup_channel_size = pileup_channel_size - 10 if no_indel else pileup_channel_size
pileup_channel_size = pileup_channel_size - 8 if no_mq else pileup_channel_size
pileup_channel_size = pileup_channel_size - 8 if no_bq else pileup_channel_size

tumor_channel_size = pileup_channel_size


flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
input_shape = [matrix_depth_dict['hifi'], no_of_positions, channel_size]
ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, channel_size]
# label_shape = [21, 3, no_of_positions, no_of_positions]


use_alt_base = True
add_af_in_label = False


label_shape = [3] if not add_af_in_label else [6]
label_size = sum(label_shape)
apply_focal_loss = True
discard_germline = False
add_l2_regulation_loss = True
somatic_arg_index = 1 if discard_germline else 2
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

min_mq = 20
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