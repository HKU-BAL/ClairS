
# ClairS

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng  

Email: rbluo@cs.hku.hk, zxzheng@cs.hku.hk  

------

## Introduction

ClairS is a deep learning-based long-read somatic variant caller.  ClairS identifies somatic SNVs (single nucleotide variants) from aligned tumor and matched normal reads efficiently and sensitively. ClairS was trained in ~30M synthetic germline datasets with various coverages and allele frequencies.

------

## Contents

- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
  - [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
- [Quick Demo](#quick-demo)
- [Usage](#usage)
- [Pre-trained Models](#pre-trained-models)
- [Disclaimer](#disclaimer)

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data, see [ONT Quick Demo](docs/ont_quick_demo.md).
- Illumina NGS data, see [Illumina Quick Demo](docs/illumina_quick_demo.md).


### General Usage

After following [installation](#installation), you can run ClairS with one command:

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont
## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

### Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available [here](https://hub.docker.com/r/hkubal/clairs). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs:latest \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## change your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \    ## change your normal bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## change your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont,ilmn}
  --output ${OUTPUT_DIR}                       ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

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
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## change your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \    ## change your normal bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## change your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont,ilmn}
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
tar -zxvf clairs_models.tar.gz -C ${CONDA_PREFIX}/bin/somatic_models/

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

## Usage

### General Usage

```bash
./run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \    ## change your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \  ## change your bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \             ## change your reference file name here
  --threads ${THREADS} \                     ## maximum threads to be used
  --platform ${PLATFORM} \                   ## options: {ont,ilmn}
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
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options {ont,ilmn}.
```

**Other parameters:**

```bash
  -P PILEUP_MODEL_PATH, --pileup_model_path PILEUP_MODEL_PATH
                        Pileup somatic model path.
  -F FULL_ALIGNMENT_MODEL_PATH, --full_alignment_model_path FULL_ALIGNMENT_MODEL_PATH
                        Full-alignment somatic model path.
  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs.
  -r REGION, --region REGION
                        Region in `ctg_name:start-end` format(1-index), default is None.
  -b BED_FN, --bed_fn BED_FN
                        Path to BED file. Call variants only in the provided bed regions.
  -V VCF_FN, --vcf_fn VCF_FN
                        Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
  -q QUAL, --qual QUAL  If set, variants with >$qual will be marked PASS, or LowQual otherwise.
  --snv_min_af SNV_MIN_AF
                        Minimum SNV AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: 0.05.
  --min_coverage MIN_COVERAGE
                        Minimum coverage required to call a somatic mutation, default: 4.
  --chunk_size CHUNK_SIZE
                        The size of each chuck for parallel processing, default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file, default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename, default: output.
  --remove_intermediate_dir
                        Remove intermediate directory, default: disable.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22} and {1..22}, default: disable.
  --print_ref_calls     Show reference calls (0/0) in VCF file, default: disable.
  --print_germline_calls
                        Show germline calls in VCF file, default: disable.
  -d, --dry_run         If true then only the commands will be printed, default: disable.
  --python PYTHON       Path of python, python3 >= 3.9 is required.
  --pypy PYPY           Path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Path of parallel, parallel >= 20191122 is required.
  --disable_phasing     If true then call variants without longphase or whatshap phasing for long-read data.
```

#### Call mutations in one or mutiple chromosomes using `-C/--ctg_name` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont -C chr21,chr22
```

#### Call mutations in one specific region using `-r/--region` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont -r chr20:1000000-2000000
```

#### Call mutations at known variant sites using `-V/--vcf_fn` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont -V input.vcf
```

#### Call mutations at multiple specific sites or bed regions using `-B/--bed_fn` parameter

We highly recommended using BED file to define multiple regions of interest like:

```shell
echo -e "${CTG1}\t${START_POS_1}\t${END_POS_1}" > input.bed
echo -e "${CTG2}\t${START_POS_2}\t${END_POS_2}" >> input.bed
...
```

Then run the command like:

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont -B input.bed
```

## Pre-trained Models

ClairS trained both pileup and full-alignment models using GIAB samples, and carry on benchmarking on HCC1395-HCC1395BL pair dataset. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

| Platform | Chemistry /Instruments | Option (`-p/--platform`) |   Reference   | Aligner  | Training samples |
| :------: | ---------------------- | ------------------------ | :-----------: | :------: | :--------------: |
|   ONT    | R10.4/R10.4.1          | `ont`                    | GRCh38_no_alt | Minimap2 |     HG001,2      |
|   ONT    | R9.4.1                 | `ont_r9`                 | GRCh38_no_alt | Minimap2 |     HG003,4      |
| Illumina | NovaSeq/HiseqX         | `ilmn`                   |    GRCh38     | BWA-MEM  |     HG003,4      |

## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.

