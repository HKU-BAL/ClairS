# Clair-somatic

## Usage

### General Usage

```bash
git clone https://github.com/zhengzhenxian/Clair-somatic.git
cd Clair-somatic

run_clair_somatic \
    --normal_bam_fn ${NORMAL_BAM_FILE_PATH} \
    --tumor_bam_fn ${TUMOR_BAM_FILE_PATH} \
    --ref_fn ${REF_FILE_PATH} \
    --output_dir ${OUTPUT_DIR} \
    --threads ${THREADS} \
    --platform ${PLATFORM}
    
## Final output file: ${OUTPUT_DIR}/somatic_output.vcf.gz
```

### 

