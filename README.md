# ClairS (early access stage)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng  
Email: rbluo@cs.hku.hk, zxzheng@cs.hku.hk  

------

## Introduction

ClairS is a somatic variant caller designed for paired samples and primarily ONT long-read. It uses Clair3 to eliminate germline variants. It ensembles the pileup and full-alignment models in Clair3, trusts them equally, and decides on the result using a set of rules and post-processing filters. With 70-fold HCC1395 (tumor) and 45-fold HCC1395BL (normal) of ONT R10.4.1 data, benchmarking against the truth SNVs (Fang et al., 2021) has shown ClairS achieved 91.10% precision and 71.92% recall. Specifically, ClairS achieved 98.12% precision and 99.15% recall on variants with AF≥0.20, and 94.26% precision and 87.46% recall on variants with AF≥0.05. Detailed performance figures are shown below.

ClairS means Clair-Somatic, or the masculine plural of "Clair" in french (thus, 's' is silent).

ClairS is now available for early access to interested and experienced users. Your suggestions and comments are highly appreciated.

------

## Performance figures

### Setup
* 70-fold HCC1395 (tumor) and 45-fold HCC1395BL (normal) of ONT R10.4.1 data.
* Benchmarking against the truths (Fang et al., 2021)
### SNV performance
<p align="center">
<img src="./images/ont_result.png" alt="ONT benmarking result" width="1080p"></img>
</p>


------

## Contents

- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
  - [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
- [Quick Demo](#quick-demo)
- [Pre-trained Models](#pre-trained-models)
- [Usage](#usage)
- [Disclaimer](#disclaimer)

------

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- Illumina NGS data as input, see [Illumina Quick Demo](docs/illumina_quick_demo.md).

### Quick start

After following [installation](#installation), you can run ClairS with one command:

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10
## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

------

## Pre-trained Models

ClairS trained both pileup and full-alignment models using GIAB samples, and carry on benchmarking on HCC1395-HCC1395BL pair dataset. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

| Platform |       Model name       | Chemistry /Instruments | Option (`-p/--platform`) |   Reference   | Aligner  | Included in the docker image | Training samples | Caveats |
| :------: | :--------------------: | :--------------------: | :----------------------: | :-----------: | :------: | :--------------------------: | :--------------: | :-----: |
|   ONT    | ont_r104_e81_sup_g5015 |     R10.4/R10.4.1      |        `ont_r10`         | GRCh38_no_alt | Minimap2 |             Yes              |     HG001,2      |         |
|   ONT    |  r941_prom_sup_g5014   |         R9.4.1         |         `ont_r9`         | GRCh38_no_alt | Minimap2 |             Yes              |     HG003,4      |    1    |
| Illumina |          ilmn          |     NovaSeq/HiseqX     |          `ilmn`          |    GRCh38     | BWA-MEM  |             Yes              |     HG003,4      |         |

**Caveats 1**: Although the r9(`r941_prom_sup_g5014`) model was trained on synthetic samples with r9.4.1 real data, the minimal AF cutoff, minimal coverage, and post-calling filtering parameters for the r9 model are copied from the r10 model, and are not optimized due to lack of real r9 data on a cancer sample with known truths.


------


## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/clairs). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs:latest \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \    ## use your normal bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10, ont_r9, ilmn}
  --output ${OUTPUT_DIR}                       ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 

```bash
conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clairs

# run the sandbox like this afterward
singularity exec clairs.sif \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs:latest \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \    ## use your normal bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10, ont_r9, ilmn}
  --output ${OUTPUT_DIR}                       ## output path prefix 
```

### Option 3. Build an anaconda virtual environment

Check here to install the tools step by step.

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install ClairS using anaconda step by step:**

```bash
# create and activate an environment named clairs
# install pypy and packages in the environemnt
conda create -n clairs -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3-illumina python=3.9.0 -y
source activate clairs

git clone https://github.com/HKU-BAL/ClairS.git
cd ClairS

# make sure in conda environment
# download pre-trained models
echo ${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin/somatic_models
wget http://www.bio8.cs.hku.hk/clairs/models/clairs_models.tar.gz
tar -zxvf clairs_models.tar.gz -C ${CONDA_PREFIX}/bin/clairs_models/

./run_clairs --help
```

### Option 4. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
git clone https://github.com/HKU-BAL/ClairS.git
cd ClairS

# build a docker image named hkubal/clairs:latest
# might require docker authentication to build docker image
docker build -f ./Dockerfile -t hkubal/clairs:latest .

# run the docker image like option 1
docker run -it hkubal/clairs:latest /opt/bin/run_clairs --help
```

------

## Usage

### General Usage

```bash
./run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \    ## use your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \  ## use your bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \             ## use your reference file name here
  --threads ${THREADS} \                     ## maximum threads to be used
  --platform ${PLATFORM} \                   ## options: {ont_r10, ont_r9, ilmn}
  --output ${OUTPUT_DIR}                     ## output path prefix
 
## Final output file: ${OUTPUT_DIR}/output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -T, --tumor_bam_fn TUMOR_BAM_FN   Tumor BAM file input. The input file must be samtools indexed.
  -N, --normal_bam_fn NORMAL_BAM_FN Normal BAM file input. The input file must be samtools indexed.
  -R, --ref_fn FASTA                Reference file input. The input file must be samtools indexed.
  -o, --output_dir OUTPUT_DIR       VCF output directory.
  -t, --threads THREADS             Max #threads to be used.
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options {ont_r10, ont_r9, ilmn}.
```

**Miscellaneous parameters:**

```bash
  -P PILEUP_MODEL_PATH, --pileup_model_path PILEUP_MODEL_PATH
                        Specify the path to your own somatic calling pileup model.
  -F FULL_ALIGNMENT_MODEL_PATH, --full_alignment_model_path FULL_ALIGNMENT_MODEL_PATH
                        Specify the path to your own somatic calling full-alignment model.
  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs. Default: all contigs will be processed.
  -r REGION, --region REGION
                        A region to be processed. Format: `ctg_name:start-end` (start is 1-based).
  -b BED_FN, --bed_fn BED_FN
                        Path to a BED file. Call variants only in the provided BED regions.
  -V VCF_FN, --vcf_fn VCF_FN
                        VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided.
  -q QUAL, --qual QUAL  If set, variants with >QUAL will be marked as PASS, or LowQual otherwise.
  --snv_min_af SNV_MIN_AF
                        Minimal SNV AF required for a variant to be called. Decrease SNV_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed and accuracy. Default: 0.05.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  --chunk_size CHUNK_SIZE
                        The size of each chuck for parallel processing. Default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename. Default: output.
  --remove_intermediate_dir
                        Remove intermediate directory before finishing to save disk space.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22} and {1..22}.
  --print_ref_calls     Show reference calls (0/0) in VCF file.
  --print_germline_calls
                        Show germline calls in VCF file.
  -d, --dry_run         Print the commands that will be ran.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --disable_phasing     Disable phasing with longphase or whatshap. Usually leads to significant performance loss.
```

#### Call SNVs in one or mutiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10 -C chr21,chr22
```

#### Call SNVs in one specific region using the `-r/--region` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10 -r chr20:1000000-2000000
```

#### Call SNVs at interested variant sites (genotyping) using the `-V/--vcf_fn` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10 -V input.vcf
```

#### Call SNVs in the BED regions using the `-B/--bed_fn` parameter

We highly recommended using BED file to define multiple regions of interest like:

```shell
echo -e "${CTG1}\t${START_POS_1}\t${END_POS_1}" > input.bed
echo -e "${CTG2}\t${START_POS_2}\t${END_POS_2}" >> input.bed
...
```

Then:

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10 -B input.bed
```

------


## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
