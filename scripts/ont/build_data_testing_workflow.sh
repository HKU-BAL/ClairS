#!/usr/bin/env bash
source /mnt/bal36/zxzheng/somatic/Clair-somatic/scripts/ont/param.sh
ROOT_FOLDER="/mnt/bal36/zxzheng/somatic/${OUTPUT_PREFIX}/test_chr20-22-workflow"
TRAIN_NAME="test"
GEN_BAMS=0

CHR=("${TEST_CHR[@]}")

TEST_FOLDER=${ROOT_FOLDER}/build
PREDICT_FOLDER=${ROOT_FOLDER}/predict
MIX_CANDIDATES_FOLDER="${TEST_FOLDER}/mixed_bam_candidates"
BAM_OUTPUT_DIR="${ROOT_FOLDER}/downsample"
BAM_INPUT_DIR="${BAM_OUTPUT_DIR}/pair"
mkdir -p ${MIX_CANDIDATES_FOLDER}
mkdir -p ${BAM_INPUT_DIR}
mkdir -p ${PREDICT_FOLDER}

if [ $GEN_CANDIDATES -eq 1 ]
then
    NORMAL_BAM_FILE_PATH=${ALL_BAM_FILE_PATH[0]}
    TUMOR_BAM_FILE_PATH=${ALL_BAM_FILE_PATH[1]}
    LOG_PATH=${BAM_OUTPUT_DIR}/logs
    DEPTH_PATH=${BAM_OUTPUT_DIR}/depth
    BAM_PATH=${BAM_OUTPUT_DIR}/split_bam
    mkdir -p ${BAM_OUTPUT_DIR}/split_bins
#    mkdir -p ${BAM_OUTPUT_DIR}/pair
    mkdir -p ${BAM_PATH}
    mkdir -p ${LOG_PATH}
    mkdir -p ${DEPTH_PATH}
#    TEST_ONLY=0
#    if [ ${TEST_ONLY} -eq 1 ]; then unset CHR; CHR=("${CHR[@]}" "${TEST_CHR[@]}"); fi
#    if [ ${TEST_ONLY} -eq 2 ]; then unset CHR; CHR=${TEST_CHR[@]}; fi
    echo ${CHR[@]}
    THREADS_LOW=20
    cd ${BAM_OUTPUT_DIR}
    #normal
    ${PARALLEL} -j${PARALLEL_THREADS} --joblog ${LOG_PATH}/parallel_1_split_chr_normal.log ${SAMTOOLS} view -@${THREADS_LOW} -bh ${NORMAL_BAM_FILE_PATH} {1} -o ${BAM_PATH}/normal_{1}.bam ::: ${CHR[@]} |& tee ${LOG_PATH}/1_SCN.log
    ${PARALLEL} -j${PARALLEL_THREADS} ${SAMTOOLS} index -@${THREADS_LOW} ${BAM_PATH}/normal_{1}.bam ::: ${CHR[@]}

    # tumor
    ${PARALLEL} -j${PARALLEL_THREADS} --joblog ${LOG_PATH}/parallel_2_split_chr_tumor.log ${SAMTOOLS} view -@${THREADS_LOW} -bh ${TUMOR_BAM_FILE_PATH} {1} -o ${BAM_PATH}/tumor_{1}.bam ::: ${CHR[@]} |& tee ${LOG_PATH}/2_SCT.log
    ${PARALLEL} -j${PARALLEL_THREADS} ${SAMTOOLS} index -@${THREADS_LOW} ${BAM_PATH}/tumor_{1}.bam ::: ${CHR[@]}
#
#    # mosdepth get depth
#    ${PARALLEL} -j${PARALLEL_THREADS} ${MOSDEPTH} -t ${THREADS} -n -x --quantize 0:15:150: ${DEPTH_PATH}/{2}_{1}.depth ${BAM_PATH}/{2}_{1}.bam ::: ${CHR[@]} ::: ${PREFIX[@]}

    # split bam
    ${PARALLEL} -j${PARALLEL_THREADS} --joblog ${LOG_PATH}/parallel_3_split_bam.log \
    ${PYPY} ${CS} SplitBam \
        --normal_bam_fn ${BAM_PATH}/normal_{1}.bam \
        --tumor_bam_fn ${BAM_PATH}/tumor_{1}.bam  \
        --output_dir ${BAM_OUTPUT_DIR}/split_bins  \
        --ctgName {1} \
        --samtools ${SAMTOOLS} \
        --cov_dir ${DEPTH_PATH} ::: ${CHR[@]} |& tee ${LOG_PATH}/3_SB.log

    # mix chunked bin
    ${PARALLEL} -j${PARALLEL_THREADS} --joblog ${LOG_PATH}/parallel_4_mix_bin.log \
    ${PYPY} ${CS} MixBin \
        --input_dir ${BAM_OUTPUT_DIR}/split_bins \
        --output_fn ${BAM_OUTPUT_DIR}/pair/tumor_{1}_{2}.bam \
        --ctgName {1} \
        --synthetic_proportion {2} \
        --synthetic_depth ${synthetic_depth} \
        --normal_bam_fn ${BAM_PATH}/normal_{1}.bam \
        --tensor_sample_mode 1 \
        --samtools ${SAMTOOLS} \
        --cov_dir ${DEPTH_PATH}  ::: ${CHR[@]} ::: ${POR[@]} |& tee ${LOG_PATH}/4_MB.log

    #softlink if no contam

    # optional mosdepth get depth
    ${PARALLEL} -j${PARALLEL_THREADS} ${MOSDEPTH} -t ${THREADS} -n -x --quantize 0:15:150: ${DEPTH_PATH}/{3}_{2}_{1}.depth ${BAM_OUTPUT_DIR}/pair/{3}_{1}_{2}.bam ::: ${CHR[@]} ::: ${POR[@]} ::: ${PREFIX[@]}
    # raw depth
    ${PARALLEL} --keep-order "echo ${DEPTH_PATH}/{2}_{1}.depth.mosdepth.summary.txt && tail -n1 ${DEPTH_PATH}/{2}_{1}.depth.mosdepth.summary.txt | cut -f 4" ::: ${CHR[@]} ::: ${PREFIX[@]} |& tee ${DEPTH_PATH}/raw_depth.log
    # sampled depth
    ${PARALLEL} --keep-order "echo ${DEPTH_PATH}/{3}_{2}_{1}.depth.mosdepth.summary.txt && tail -n1 ${DEPTH_PATH}/{3}_{2}_{1}.depth.mosdepth.summary.txt | cut -f 4" ::: ${CHR[@]} ::: ${POR[@]} ::: ${PREFIX[@]} |& tee ${DEPTH_PATH}/sampled_depth.log

fi

mkdir -p ${MIX_CANDIDATES_FOLDER}
cd ${MIX_CANDIDATES_FOLDER}
NORMAL_ALT_PATH="${MIX_CANDIDATES_FOLDER}/normal_alt_output"
TUMOR_ALT_PATH="${MIX_CANDIDATES_FOLDER}/tumor_alt_output"
LOG_PATH="${MIX_CANDIDATES_FOLDER}/log"
PLOT_PATH="${MIX_CANDIDATES_FOLDER}/plot"
REF_CANS_PATH="${MIX_CANDIDATES_FOLDER}/ref_cans"
mkdir -p ${MIX_CANDIDATES_FOLDER}/tmp
mkdir ${NORMAL_ALT_PATH}
mkdir ${TUMOR_ALT_PATH}
mkdir ${LOG_PATH}
mkdir ${PLOT_PATH}
mkdir ${REF_CANS_PATH}


CANDIDATES_FOLDER_PATH="${TEST_FOLDER}/candidates"
SPLIT_BED_PATH="${TEST_FOLDER}/split_beds"
LOG_PATH="${TEST_FOLDER}/log"
NORMAL_TENSOR_PATH="${TEST_FOLDER}/normal_tensor_output"
NORMAL_ALT_PATH="${TEST_FOLDER}/normal_alt_output"
TUMOR_TENSOR_PATH="${TEST_FOLDER}/tumor_tensor_output"
BINS_FOLDER_PATH="${TEST_FOLDER}/bins"
UNIFIED_VCF_PATH="${TEST_FOLDER}/unified_vcf"
mkdir -p ${TEST_FOLDER}
mkdir -p ${TEST_FOLDER}/tmp
mkdir -p ${TEST_FOLDER}/tensor_can
mkdir -p ${PREDICT_FOLDER}/predict
mkdir -p ${PREDICT_FOLDER}/vcf_output
mkdir -p ${CANDIDATES_FOLDER_PATH}
#mkdir -p ${CANDIDATES_FOLDER_PATH}/fp
#mkdir -p ${CANDIDATES_FOLDER_PATH}/tp
mkdir -p ${NORMAL_TENSOR_PATH}
mkdir -p ${TUMOR_TENSOR_PATH}
mkdir -p ${NORMAL_ALT_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${UNIFIED_VCF_PATH}
mkdir -p ${LOG_PATH}

MODEL_FOLDER_PATH="${ROOT_FOLDER}/train/${TRAIN_NAME}"
mkdir -p ${MODEL_FOLDER_PATH}/validate


cd ${MIX_CANDIDATES_FOLDER}

    echo "[INFO] Extract Normal Candidates"
    time ${PARALLEL} --joblog ${LOG_PATH}/parallel_1_extract_candidates.log -j${THREADS} \
    ${PYPY} ${CS} ExtractAF \
        --bam_fn ${BAM_INPUT_DIR}/${PREFIX[0]}_{1/.}_{7}.bam \
        --ref_fn {5} \
        --ctgName {1} \
        --samtools ${SAMTOOLS} \
        --min_af ${MIN_AF} \
        --chunk_id {8} \
        --chunk_num ${chunk_num} \
        --alt_fn ${NORMAL_ALT_PATH}/{2}_{3}_{1}_{7}_{8} \
        --test_pos False \
        --output_alt_info \
        --extend_bed ${INTERSECTED_BED} \
        --output_depth ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BAM_FILE_PATH[0]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} :::+ ${PRO[@]} ::: ${CHUNK_LIST[@]} |& tee ${LOG_PATH}/1_EC.log

echo "[INFO] Extract Tumor Candidates"
time ${PARALLEL} --joblog ${LOG_PATH}/parallel_2_extract_tumor_candidates.log -j${THREADS} \
${PYPY} ${CS} extract_candidates \
    --bam_fn ${TUMOR_BAM_FILE_PATH} \
    --ref_fn {5} \
    --ctgName {1} \
    --samtools ${SAMTOOLS} \
    --snp_min_af ${MIN_SNP_AF} \
    --indel_min_af ${MIN_INDEL_AF} \
    --vcf_fn ${UNIFIED_VCF[1]} \
    --chunk_id {8} \
    --chunk_num ${chunk_num} \
    --alt_fn ${TUMOR_ALT_PATH}/{2}_{3}_{1}_{7}_{8} \
    --test_pos False \
    --output_alt_info \
    --bed_fn ${INTERSECTED_BED} \
    --extend_bed ${INTERSECTED_BED} \
    --candidates_folder ${CANDIDATES_FOLDER_PATH} \
    --output_depth ::: ${CHR[@]} ::: ${ALL_SAMPLE[1]} :::+ ${DEPTHS[1]} :::+ ${ALL_BAM_FILE_PATH[1]} :::+ ${ALL_REFERENCE_FILE_PATH[1]} :::+ ${ALL_BED_FILE_PATH[1]} ::: ${POR[@]} ::: ${CHUNK_LIST[@]} |& tee ${LOG_PATH}/2_ETC.log
#   --min_truth_af ${MIN_TRUTH_AF} \
#   --store_tumor_infos \

${PARALLEL} -j${THREADS} "gzip -fdc {2} | grep -v '#' | grep -w {1} > ${UNIFIED_VCF_PATH}/unified_{3}_{1}" ::: ${CHR[@]} ::: ${UNIFIED_VCF[@]} :::+ ${ALL_SAMPLE[@]}

echo "[INFO] Get Candidates"
${PARALLEL} -j${THREADS} --joblog ${LOG_PATH}/parallel_1_get_candidates.log \
${PYPY} ${CS}  GetCandidates \
    --vcf_fn_1 ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[0]}_{1} \
    --vcf_fn_2 ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[1]}_{1}  \
    --normal_reference_cans ${MIX_CANDIDATES_FOLDER}/normal_alt_output/${ALL_SAMPLE[0]}_${DEPTHS[0]}_{1}_ \
    --tumor_reference_cans ${MIX_CANDIDATES_FOLDER}/tumor_alt_output/${ALL_SAMPLE[1]}_${DEPTHS[1]}_{1}_ \
    --bed_fn ${INTERSECTED_BED}  \
    --ctgName {1} \
    --split_folder ${CANDIDATES_FOLDER_PATH} ::: ${CHR[@]} |& tee ${LOG_PATH}/1_GC.log

cat ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILE_* > ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILES

echo "[INFO] Extract Pair Candidates"
time ${PARALLEL}  --joblog ${LOG_PATH}/parallel_2_create_pair_tensor.log -j${THREADS} \
${PYPY} ${CS} create_pair_tensor \
    --normal_bam_fn ${NORMAL_BAM_FILE_PATH} \
    --tumor_bam_fn ${TUMOR_BAM_FILE_PATH} \
    --ref_fn ${ALL_REFERENCE_FILE_PATH[0]} \
    --ctgName {1/.} \
    --samtools ${SAMTOOLS} \
    --full_aln_regions {1} \
    --tensor_can_fn ${TEST_FOLDER}/tensor_can/{2}_{1/} \
    --alt_fn ${NORMAL_ALT_PATH}/{2}_{1/} \
    --platform ${platform} \
    --truth_vcf_fn ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[1]}_{1/.} \
    --test_pos False :::: ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILES ::: ${POR[@]} |& tee ${LOG_PATH}/2_CPT.log

GPU_ID=`nvidia-smi | tail -n 2 | head -n 1`
if [ "$GPU_ID" == '|  No running processes found                                                 |' ]; then GPU_ID="0"; else GPU_ID="1"; fi
echo "Avaliable GPU ID:${GPU_ID}"
#
#if [ ${GPU_SERVER["$HOSTNAME"_0]} = '3090' ]
#then
#    if [ "$GPU_ID" == "0" ]
#    then
#        PYTHON_TORCH="/autofs/bal33/zxzheng/env/miniconda2/envs/torch_3090/bin/python3"
#    else
#        PYTHON_TORCH="/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3"
#    fi
#else
#    if [ "$GPU_ID" == "0" ]
#    then
#        PYTHON_TORCH="/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3"
#    else
#        echo -e "${ERROR} NO GPU avaliable${NC}"
#        PYTHON_TORCH="/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3"
#    fi
#fi
#echo "PYTHONPATH: ${PYTHON_TORCH}"

use_gpu=0
if [ ${use_gpu} -eq 1 ]
then
PREDICT_THREADS=6
else
PREDICT_THREADS=36
fi

time ${PARALLEL} --joblog ${TEST_FOLDER}/parallel_predict_log.txt -j${PREDICT_THREADS} ${PYTHON_TORCH} ${CS} predict \
--tensor_fn ${TEST_FOLDER}/tensor_can/{2}_{1/} \
--predict_fn ${PREDICT_FOLDER}/predict/{2}_{1/} \
--chkpnt_fn /mnt/bal36/zxzheng/somatic/ont/12345/train/torch_test/9.pkl \
--use_gpu ${use_gpu} \
--platform ${platform} :::: ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILES ::: ${POR[@]} |& tee ${TEST_FOLDER}/PREDICT.log


#--add_indel_length ${add_indel_length} \
time ${PARALLEL} --joblog ${TEST_FOLDER}/parallel_call_variants_log.txt -j${PREDICT_THREADS} ${PYTHON_TORCH} ${CS} call_variants \
--predict_fn ${PREDICT_FOLDER}/predict/{2}_{1/} \
--call_fn ${PREDICT_FOLDER}/vcf_output/{2}_{1/}.vcf \
--ref_fn ${ALL_REFERENCE_FILE_PATH[0]} \
--platform ${platform} :::: ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILES ::: ${POR[@]} |& tee ${TEST_FOLDER}/CV.log

cat ${TEST_FOLDER}/vcf_output/*.vcf | ${PYPY} ${CS} sort_vcf --ref_fn ${ALL_REFERENCE_FILE_PATH[0]} --output_fn ${TEST_FOLDER}/output.vcf

${PYPY} ${CLAIR3} SortVcf \
    --input_dir ${FULL_ALIGNMENT_OUTPUT_PATH} \
    --vcf_fn_prefix "full_alignment" \
    --output_fn ${OUTPUT_FOLDER}/full_alignment.vcf \
    --sampleName ${SAMPLE} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --contigs_fn ${TMP_FILE_PATH}/CONTIGS

#--sampleName {2} \
#--chunk_id {6} \
#--is_from_tables \
#--chunk_num ${call_chunk_num} \
#--predict_fn ${OUTPUT_ROOT_DIRECTORY}/${TRAIN_NAME}/epoch{4}/predict_{2}_{3}_{1}_{5}_{6} \


USE_RESNET=0
ADD_WRITER=0
MODEL_FOLDER_PATH=${MODEL_FOLDER_PATH}_parallel
mkdir -p ${MODEL_FOLDER_PATH}/validate
cd ${MODEL_FOLDER_PATH}

PREDICT_THREADS=3
call_chunk_num=3
CALL_CHUNK_LIST=`seq 1 ${call_chunk_num}`

export CUDA_VISIBLE_DEVICES="${GPU_ID}"
echo "[INFO] Model testing"
time ${PARALLEL} --joblog ${MODEL_FOLDER_PATH}/parallel_predict.log -j${PREDICT_THREADS} \
${PYTHON_TORCH} ${CS} Predict \
    --bin_fn ${BINS_FOLDER_PATH} \
    --platform ${platform} \
    --chkpnt_fn /mnt/bal36/zxzheng/somatic/ont/12345/train/torch_test/18.pkl \
    --output_logits \
    --output_dir ${MODEL_FOLDER_PATH}/validate \
    --ctgName ${TEST_CHR_S} \
    --test_all_pos \
    --use_resnet ${USE_RESNET} \
    --add_writer ${ADD_WRITER} \
    --unified_vcf_fn ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[1]}_ \
    --chunk_id {1} \
    --chunk_num ${call_chunk_num} ::: ${CALL_CHUNK_LIST[@]} |& tee ${MODEL_FOLDER_PATH}/PREDICT.log


PYTHON_TORCH="/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3"
time ${PARALLEL} --joblog ${PREDICT_FOLDER}/parallel_predict.log -j${THREADS} \
${PYTHON_TORCH} ${CS} call_variants_from_bam \
    --normal_bam_fn ${NORMAL_BAM_FILE_PATH} \
    --tumor_bam_fn ${TUMOR_BAM_FILE_PATH} \
    --ref_fn ${ALL_REFERENCE_FILE_PATH[0]} \
    --ctgName {1/.} \
    --samtools ${SAMTOOLS} \
    --full_aln_regions {1} \
    --platform ${platform} \
    --python /mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/python3 \
    --chkpnt_fn /mnt/bal36/zxzheng/somatic/ont/12345/train/torch_test/9.pkl \
    --call_fn ${PREDICT_FOLDER}/vcf_output/{2}_{1/}.vcf :::: ${CANDIDATES_FOLDER_PATH}/FULL_ALN_FILES ::: ${POR[@]} |& tee ${PREDICT_FOLDER}/PREDICT.log
#    --truth_vcf_fn ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[1]}_{1/.} \



#
#chunk_num=30
#CHUNK_LIST=`seq 1 ${chunk_num}`
#PYPY="/mnt/bal36/zxzheng/env/miniconda3/envs/clair3/bin/pypy3"
#FP_FN_LIST=(fp fn tp)
#mkdir -p ${MODEL_FOLDER_PATH}/fp/backup
#mkdir -p ${MODEL_FOLDER_PATH}/fn/backup
#mkdir -p ${MODEL_FOLDER_PATH}/tp/backup
#${PARALLEL} -j${THREADS} \
#${PYPY} ${CS} PlotAlignment \
#    --normal_bam_fn ${BAM_INPUT_DIR}/${PREFIX[0]}_{1}_{2}.bam \
#    --tumor_bam_fn ${BAM_INPUT_DIR}/${PREFIX[1]}_{1}_{2}.bam \
#    --raw_tumor_bam_fn ${ALL_BAM_FILE_PATH[1]} \
#    --raw_normal_bam_fn ${ALL_BAM_FILE_PATH[0]} \
#    --ref_fn ${ALL_REFERENCE_FILE_PATH[0]} \
#    --ctgName {1} \
#    --sampleName 'hg002' \
#    --output_dir ${MODEL_FOLDER_PATH}/{4} \
#    --unified_vcf_fn ${MODEL_FOLDER_PATH}/validate/{4}.vcf \
#    --truth_vcf_fn ${UNIFIED_VCF_PATH}/unified_${ALL_SAMPLE[1]}_{1} \
#    --chunk_id {3} \
#    --test_pos False \
#    --bed_fn ${INTERSECTED_BED} \
#    --chunk_num ${chunk_num} ::: ${CHR[@]} ::: ${POR[@]} ::: ${CHUNK_LIST[@]} ::: ${FP_FN_LIST[@]}

#    --use_siam \
#    --add_contrastive
# --exclude_training_samples ${exclude_training_samples}
#--add_indel_length ${add_indel_length} \
#/autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/python3 /mnt/bal36/zxzheng/somatic/Clair-somatic/clair-somatic.py Predict --bin_fn /mnt/bal36/zxzheng/somatic/ont/test_chr20/build/bins --platform ont --chkpnt_fn /mnt/bal36/zxzheng/somatic/ont/all/train/test/14 --output_logits --output_dir /mnt/bal36/zxzheng/somatic/ont/test_chr20/train/test/validate --ctgName chr20 --test_all_pos --unified_vcf_fn /mnt/bal36/zxzheng/somatic/ont/test_chr20/build/unified_vcf/unified_hg004_