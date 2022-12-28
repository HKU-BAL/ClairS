# Parameters
INPUT_DIR="${HOME}/illumina_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38 reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/GRCh38_chr17.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/GRCh38_chr17.fa.fai
# Normal and tumor BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395BL_normal_chr17_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395BL_normal_chr17_demo.bam.bai
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395_tumor_chr17_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395_tumor_chr17_demo.bam.bai

# SEQC2 Truth VCF and BED
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_High-Confidence_Regions_v1.2_chr17.bed

REF="GRCh38_chr17.fa"
NORMAL_BAM="HCC1395BL_normal_chr17_demo.bam"
TUMOR_BAM="HCC1395_tumor_chr17_demo.bam"
BASELINE_VCF_FILE_PATH="SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz"
BASELINE_BED_FILE_PATH="SEQC2_High-Confidence_Regions_v1.2_chr17.bed"
OUTPUT_VCF_FILE_PATH="output.vcf.gz"

# Run clair-somatic using one command
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-somatic:latest \
  /opt/bin/run_clair_somatic \
  --tumor_bam_fn ${INPUT_DIR}/${TUMOR_BAM} \
  --normal_bam_fn ${INPUT_DIR}/${NORMAL_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform ilmn \
  --output_dir ${OUTPUT_DIR} \
  --region chr17:80000000-80100000

#Benchmarking
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-somatic:latest \
  python3 /opt/bin/clair-somatic.py compare_vcf \
     --truth_vcf_fn ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
     --input_vcf_fn ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
     --bed_fn ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
     --output_dir ${OUTPUT_DIR}/benchmark \
     --input_filter_tag 'PASS' \
     --ctg_name chr17 \
     --ctg_start 80000000 \
     --ctg_end 80100000
