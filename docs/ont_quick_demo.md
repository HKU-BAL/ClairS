## ONT Somatic Calling Quick Demo
Here is a quick demo for the Oxford Nanopore (ONT) somatic calling using HCC1395 tumor-normal pair chromosome 17 data. The data was sequenced using ONT [Q20+](https://nanoporetech.com/q20plus-chemistry) chemistry.

```bash
Platform:          ONT
Sample:     	   HCC1395-HCC1395BL tumor-normal pair
Normal coverage:   ~40x
Tumor coverage:    ~70x
Reference:         GRCh38_no_alt
Aligner:           minimap2
Region:            chr17:80000000-80100000
Basecaller:        Guppy 6.1.5
Chemistry:         R10.4.1
```

**Download data**

```bash
# Parameters
INPUT_DIR="${HOME}/ont_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/GRCh38_no_alt_chr17.fa
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/GRCh38_no_alt_chr17.fa.fai
# Normal and tumor BAM
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/HCC1395BL_normal_chr17_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/HCC1395BL_normal_chr17_demo.bam.bai
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/HCC1395_tumor_chr17_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ont/HCC1395_tumor_chr17_demo.bam.bai

# SEQC2 Truth VCF and BED
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair_somatic/quick_demo/ilmn/SEQC2_High-Confidence_Regions_v1.2_chr17.bed

REF="GRCh38_no_alt_chr17.fa"
NORMAL_BAM="HCC1395BL_normal_chr17_demo.bam"
TUMOR_BAM="HCC1395_tumor_chr17_demo.bam"
BASELINE_VCF_FILE_PATH="SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz"
BASELINE_BED_FILE_PATH="SEQC2_High-Confidence_Regions_v1.2_chr17.bed"
OUTPUT_VCF_FILE_PATH="output.vcf.gz"

```

#### Somatic variant calling using docker pre-built image

```bash
# Run clair-somatic using one command
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-somatic:latest \
  /opt/bin/run_clair_somatic.sh \
  --tumor_bam_fn ${INPUT_DIR}/${TUMOR_BAM} \
  --normal_bam_fn ${INPUT_DIR}/${NORMAL_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform ont \
  --output ${OUTPUT_DIR} \
  --region chr17
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
| **SNV** |    1.0    | 0.931  |  0.9643  |  27  |  0   |  2   |

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
wget "https://raw.githubusercontent.com/HKU-BAL/Clair-Somatic/main/demo/ont_quick_demo.sh"
chmod +x ont_quick_demo.sh
./ont_quick_demo.sh
```

Check the results using `less ${HOME}/ont_quick_demo/output/output.vcf.gz`.

