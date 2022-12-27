## Illumina Somatic Variant Calling Quick Demo
Here is a quick demo for the Illumina NGS soamtic variant calling using somatic calling using HCC1395 tumor-normal pair chromosome 17 data. The data was acquired from [SEQC2](https://sites.google.com/view/seqc2/home?authuser=0).

```bash
Platform:          Illumina
Sample:            HCC1395-HCC1395BL tumor-normal pair
Normal coverage:   ~50x
Tumor coverage:    ~50x
Reference:         GRCh38
Aligner:           BWA-MEM
Region:            chr17:80000000-80100000
Instruments:       NovaSeq 6000
```

### Download data

```bash
# Parameters
INPUT_DIR="${HOME}/illumina_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38 reference
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/GRCh38_chr17.fa
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/GRCh38_chr17.fa.fai
# Normal and tumor BAM
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395BL_normal_chr17_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395BL_normal_chr17_demo.bam.bai
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395_tumor_chr17_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/HCC1395_tumor_chr17_demo.bam.bai

# SEQC2 Truth VCF and BED
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_High-Confidence_Regions_v1.2_chr17.bed

REF="GRCh38_chr17.fa"
NORMAL_BAM="HCC1395BL_normal_chr17_demo.bam"
TUMOR_BAM="HCC1395_tumor_chr17_demo.bam"
BASELINE_VCF_FILE_PATH="SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz"
BASELINE_BED_FILE_PATH="SEQC2_High-Confidence_Regions_v1.2_chr17.bed"
OUTPUT_VCF_FILE_PATH="output.vcf.gz"

```

### Somatic variant calling using docker pre-built image

```bash
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
  --output ${OUTPUT_DIR} \
  --region chr17:80000000-80100000
```

**Run [compare_vcf.py](src/compare.vcf) for benchmarking (optional)**

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-somatic:latest \
  python3 /opt/bin/src/compare_vcf.py \
     --truth_vcf_fn ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
     --input_vcf_fn ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
     --bed_fn ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
     --output_dir ${OUTPUT_DIR}/benchmark \
     --input_filter_tag 'PASS' \
     --ctg_name chr17
```

**Expected output:**

|  Type   | Precision | Recall | F1-score |  TP  |  FP  |  FN  |
| :-----: | :-------: | :----: | :------: | :--: | :--: | :--: |
| **SNV** |    1.0    |  1.0   |   1.0    |  29  |  0   |  0   |

 **Or run [som.py]() for benchmarking (optional)**

```bash
# Run hap.py
docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/som.py \
    ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
    -T ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
    -f ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
    -r ${INPUT_DIR}/${REF} \
    -o "${OUTPUT_DIR}/som" \
    -l chr17
```

**Run all commands above:**

```bash
cd ${HOME}
wget "https://raw.githubusercontent.com/HKU-BAL/Clair-Somatic/main/demo/ilmn_quick_demo.sh"
chmod +x ilmn_quick_demo.sh
./ilmn_quick_demo.sh
```

Check the results using `less ${HOME}/illumina_quick_demo/output/output.vcf.gz`.