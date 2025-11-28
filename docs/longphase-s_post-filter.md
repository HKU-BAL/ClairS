## Enhancing ClairS Somatic Calling with LongPhase-S

[LongPhase-S](https://github.com/CCU-Bioinformatics-Lab/longphase-s) is a haplotype phasing tool that  jointly reconstructs somatic haplotypes, infers tumor purity, and recalibrates somatic variants in a purity-aware manner for paired tumor-normal long-read sequencing. This workflow demonstrates how to integrate LongPhase-S with ClairS to improve somatic variant detection through haplotype-aware analysis.

The workflow consists of three steps:
1. Run ClairS to identify somatic and germline variants
2. Use LongPhase-S to phase germline variants from the normal sample
3. Apply phasing information to haplotag tumor reads and refine somatic calls

---

### Step 1: Run ClairS for Somatic Variant Calling

First, run ClairS on the tumor-normal pair to identify somatic variants and germline variants:
```
# Set input parameters
INPUT_DIR="${HOME}/input"
OUTPUT_DIR="${HOME}/output"
REF="${INPUT_DIR}/reference.fa"
NORMAL_BAM="${INPUT_DIR}/normal.bam"
TUMOR_BAM="${INPUT_DIR}/tumor.bam"
THREADS=32
PLATFORM="ont_r10_dorado_sup_5khz"

# Run ClairS
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs:v0.4.4 \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${TUMOR_BAM} \
  --normal_bam_fn ${NORMAL_BAM} \
  --ref_fn ${REF} \
  --threads ${THREADS} \
  --platform ${PLATFORM} \
  --output_dir ${OUTPUT_DIR} \
  --sample_name "sample" **Key outputs:**
```
- `${OUTPUT_DIR}/output.vcf.gz`: Somatic variants called by ClairS
- `${OUTPUT_DIR}/tmp/clair3_output/clair3_normal_output/merge_output.vcf.gz`: Germline SNPs called by Clair3 (used for phasing)

---

### Step 2: Phase Germline Variants with LongPhase-S

Phase the germline variants from the normal sample using LongPhase-S:
```
# Create output directory
mkdir -p ${OUTPUT_DIR}/longphase_s_output

GERMLINE_VCF="${OUTPUT_DIR}/tmp/clair3_output/clair3_normal_output/merge_output.vcf.gz"

longphase-s phase \
  -s ${GERMLINE_VCF} \
  -b ${NORMAL_BAM} \
  -r ${REF} \
  -t ${THREADS} \
  -o ${OUTPUT_DIR}/longphase_s_output/phased \
  --ont
```
- `${OUTPUT_DIR}/longphase_s_output/phased.vcf`: Phased germline variants with haplotype information

---

### Step 3: Haplotag Tumor BAM and Refine Somatic Calls

Use the phased germline variants to haplotag tumor reads and refine somatic variant calls:
```
longphase-s somatic_haplotag \
  -s ${OUTPUT_DIR}/longphase_s_output/phased.vcf \
  -b ${NORMAL_BAM} \
  --tumor-snp-file ${OUTPUT_DIR}/output.vcf.gz \
  --tumor-bam-file ${TUMOR_BAM} \
  -r ${REF} \
  -t ${THREADS} \
  -o ${OUTPUT_DIR}/longphase_s_output/haplotagged_tumor \
  --tagSupplementary \
  -q 20 \
  --output-somatic-vcf
```
- `${OUTPUT_DIR}/longphase_s_output/haplotagged_tumor.bam`: Haplotagged tumor BAM with HP tags
- `${OUTPUT_DIR}/longphase_s_output/haplotagged_tumor_somatic.vcf`: Refined somatic variants with phasing information

**Key parameters:**
- `--ont`: Use ONT-specific parameters for phasing
- `--tagSupplementary`: Include supplementary alignments in haplotagging
- `-q 20`: Minimum mapping quality threshold
- `--output-somatic-vcf`: Generate refined somatic VCF with haplotype information

---
