<div align="center">
    <img src="images/clairs_icon.png" width = "200" alt="ClairS">
</div>

# ClairS - a deep-learning method for long-read somatic small variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng  
Email: rbluo@cs.hku.hk, zxzheng@cs.hku.hk  

------

## Introduction

ClairS is a somatic variant caller designed for paired samples and primarily ONT long-read. It uses Clair3 to eliminate germline variants. It ensembles the pileup and full-alignment models in Clair3, trusts them equally, and decides on the result using a set of rules and post-processing filters. With 50-fold HCC1395 (tumor) and 25-fold HCC1395BL (normal) of ONT R10.4.1 data, benchmarking against the truth SNVs ([Fang et al., 2021](https://www.nature.com/articles/s41587-021-00993-6)), ClairS achieved 86.86%/93.01% recall/precision rate for SNVs when targeting VAF ≥0.05. For variants with VAF ≥0.2, the numbers go up to 94.65%/96.63%. Detailed performance figures are shown below.

ClairS means Clair-Somatic, or the masculine plural of "Clair" in french (thus, 's' is silent).

The logo of ClairS was generated using DALL-E 2 with prompt "A DNA sequence with genetic variant that looks like a letter 'S'".

A preprint describing ClairS's algorithms and results is at [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.08.17.553778v1).

For germline variant calling using **DNA-seq** sample, please try [Clair3](https://github.com/HKU-BAL/Clair3). 

For germline variant calling using **long-read RNA-seq** sample, please try [Clair3-RNA](https://github.com/HKU-BAL/Clair3-RNA).

For somatic variant calling using **tumor only** sample, please try [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO).

------

## Performance figures

### ONT Q20+ chemistry performance
The latest performance figures as of Oct 10th, 2024 (ClairS v0.4.0) is available in this [technical note](docs/Improving_the_performance_of_ClairS_and_ClairS-TO_with_new_real_cancer_cell-line_datasets_and_PoN.pdf).

Performance comparison between “ClairS v0.4.0 with the SS model”, “ClairS v0.4.0 with the SS+RS model”, and “DeepSomatic v1.7.0”, at (a) different coverages, and (b) at different AF ranges for SNV and Indel, respectively.  
  
![](./images/ont_performance.png)  

### PacBio Revio SNV performance

- HCC1395/HCC1395BL tumor/normal of PacBio Revio data, using SMRTbell prep kit 3.0
- Truth:High confidence (HighConf) and medium confidence (MedConf) SNV from the SEQC2 HCC1395/BL truths ([Fang et al., 2021](https://www.nature.com/articles/s41587-021-00993-6)), the TVAF (tumor variant allele frequency) of which is ≥0.05 in the above dataset

#### The performance of ClairS at multiple VAF ranges and multiple tumor coverages with the normal coverage fixed at 25x

![](./images/hifi_vaf_1_result.png)

#### The performance of ClairS at multiple VAF ranges and multiple normal coverages with the tumor coverage fixed at 50x

![](./images/hifi_vaf_2_result.png)



### Illumina SNV performance

- HCC1395/HCC1395BL tumor/normal of of Illumina NovaSeq 6000 and HiSeq 4000 data
- Truth:High confidence (HighConf) and medium confidence (MedConf) SNV from the SEQC2 HCC1395/BL truths ([Fang et al., 2021](https://www.nature.com/articles/s41587-021-00993-6)), the TVAF (tumor variant allele frequency) of which is ≥0.05 in the above dataset

#### The precision-recall curve of different tumor/normal purity combinations with tumor coverage fixed at 50x and normal coverage fixed at 25x

![](./images/illumina_pr_curve_result.png)



------

## Contents
- [Latest Updates](#latest-updates)
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

## Latest Updates
*v0.4.2 (Jue 29, 2025)* : Added `--snv_min_qual` and `--indel_min_qual` options to independently set the minimum QUAL threshold for SNVs and Indels to be marked as 'PASS', while deprecating the legacy `--qual` option.

*v0.4.1 (Nov 29)* : Added ssrs model for PacBio Revio (`hifi_revio_ssrs`) and illumina (`ilmn_ssrs`) platforms.

*v0.4.0 (Oct 11)* : This version is a major update. The new features and benchmarks are explained in a technical note titled [“Improving the performance of ClairS and ClairS-TO with new real cancer cell-line datasets and PoN”](docs/Improving_the_performance_of_ClairS_and_ClairS-TO_with_new_real_cancer_cell-line_datasets_and_PoN.pdf). A summary of changes: 1. Starting from this version, ClairS will provide two model types. `ssrs` is a model trained initially with synthetic samples and then real samples augmented (e.g., `ont_r10_dorado_sup_5khz_ssrs`), `ss` is a model trained from synthetic samples (e.g., `ont_r10_dorado_sup_5khz_ss`). The `ssrs` model provides better performance and fits most usage scenarios. `ss` model can be used when missing a cancer-type in model training is a concern. In v0.4.0, four real cancer cell-line datasets (HCC1937/BL, HCC1954/BL, H1437/BL, and H2009/BL) covering two cancer types (breast cancer, lung cancer) published by [Park et al.](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1) were used for `ssrs` model training. 2. Added BQ jittering in model training to address the BQ distribution difference between the training and calling datasets that leads to performance drop. 3. Added the `--indel_min_af` option and adjusted the default minimum allelic fraction requirement to 0.1 for Indels in ONT platform.

*v0.3.1 (Aug 16)* : 1. Added four options i. `--use_heterozygous_snp_in_tumor_sample_and_normal_bam_for_intermediate_phasing`, ii. `--use_heterozygous_snp_in_normal_sample_and_normal_bam_for_intermediate_phasing`, iii. `--use_heterozygous_snp_in_tumor_sample_and_tumor_bam_for_intermediate_phasing`, and iv. `--use_heterozygous_snp_in_normal_sample_and_tumor_bam_for_intermediate_phasing`. iii is equivalent to `--use_heterozygous_snp_in_tumor_sample_for_intermediate_phasing` added in v0.2.0. iv is equivalent to `--use_heterozygous_snp_in_normal_sample_for_intermediate_phasing` added in v0.2.0. Use normal bam for intermediate phasing was a request from @[Sergey Aganezov](https://github.com/aganezov). When the coverage of normal and tumor are similar, using normal bam for intermediate phasing has negligible difference from using tumor bam in our experiments using HCC1395/BL. 2. Added `--haplotagged_tumor_bam_provided_so_skip_intermediate_phasing_and_haplotagging` to use the haplotype information provided in the tumor bam directly and skip intermediate phasing and haplotagging. This option is useful when using ClairS in a pipeline in which the phasing of the tumor bam is done before running ClairS. BAM haplotagged by WhatsHap and LongPhase are accepted. 3. Bumped up Clair3 dependency to version 1.0.10, LongPhase to version 1.7.3.

*v0.3.0 (Jul 5)* : 1. Added a module called “verdict” (Option `--enable_verdict`) to statistically classify a called variant into either a germline, somatic, or subclonal somatic variant based on the CNV profile and tumor purity estimation. Please find out more technical details about the Verdict module [here](docs/verdict.md). 2. Improved model training speed, reduced model training time cost by about three times.

*v0.2.0 (Apr 29)* : 1. Added `--use_heterozygous_snp_in_normal_sample_for_intermediate_phasing`/`--use_heterozygous_snp_in_tumor_sample_for_intermediate_phasing` option to support using either heterozygous SNPs in the normal sample or tumor sample for intermediate phasing. The previous versions used in_tumor_sample for phasing. In this new version, when testing with ONT 4kkz HCC1395/BL and using in_normal_sample for intermediate phasing, the SNV precision improved ~2%, while recall remained unchanged. in_normal_sample becomes the default from this version. However, if the coverage of normal sample is low, please consider switching back to using in_tumor_sample ([#22](https://github.com/HKU-BAL/ClairS/issues/22), idea contributed by the longphase team @[sloth-eat-pudding](https://github.com/sloth-eat-pudding)). 2. Added `--use_heterozygous_indel_for_intermediate_phasing` to include high quality heterozygous Indels for intermediate phasing. With this new option, the haplotagged tumor reads increased by ~3% in ONT 4khz HCC1395/BL, the option becomes default from this version. 3. Added a model that might provide a slightly better performance for liquid tumor. In this release, only ONT Dorado 5khz HAC for liquid tumor (`-p ont_r10_dorado_hac_5khz_liquid`) is provided. The model was trained with slightly higher normal contamination. We are testing out the new model with collaborator. 4. Added `--use_longphase_for_intermediate_haplotagging` option to replace WhatsHap haplotagging by LongPhase haplotagging to speed up read haplotagging process, the option becomes default from this version. 5. Bumped up Clair3 dependency to version 1.0.7, LongPhase to version 1.7.

*v0.1.7 (Jan 25, 2024)* : 1. Added ONT Dorado 5khz HAC (`-p ont_r10_dorado_hac_5khz`) and Dorado 4khz HAC (`-p ont_r10_dorado_hac_4khz`) model, renamed all ONT Dorado SUP model, check [here](https://github.com/HKU-BAL/ClairS/blob/main/README.md#pre-trained-models) for more details. 2. Enabled somatic variant calling in sex chromosomes. 3. Added `FAU`, `FCU`, `FGU`, `FTU`, `RAU`, `RCU`, `RGU`, and `RTU` tags.

*v0.1.6 (Sep 18)* : 1. Fixed an output bug that caused no VCF output if no Indel candidate was found (contributor @[Khi Pin](https://github.com/proteinosome)). 2. Fixed showing incorrect reference allele depth at a deletion region. 3. Added PacBio HiFi [quick demo](docs/pacbio_hifi_quick_demo.md).

*v0.1.5 (Aug 2)* : 1. Updated SNV calling using ONT Dorado 4kHz data with a new model trained using multiple-sample pairs (HG003/4); 2. Updated SNV calling using ONT Dorado 5kHz data with a new model trained using multiple-sample pairs (HG001/HG002, HG003/4); 3.  Support somatic indel calling using ONT Dorado 4kHz data. 4. Support somatic indel calling using ONT Dorado 5kHz data.

*v0.1.4 (Jul 15)* : 1. Added reference depth in AD tag. 2. Added HiFi Sequel II Indel model.

*v0.1.3 (Jul 5)* : Added ONT Dorado 4khz (`-p ont_r10_dorado_4khz`) and 5khz (`-p ont_r10_dorado_5khz`) models, check [here](#pre-trained-models) for more details. Renamed platform options `ont_r10` to `ont_r10_guppy` and `ont_r9` to `ont_r9_guppy`.

*v0.1.2 (May 17)* : Added HiFi Revio model, renamed HiFi Sequel II model from `hifi` to `hifi_sequel2`.

*v0.1.1 (Apr 30)* : 1. Added the "command line used" to VCF header. 2. Added `NAU`, `NCU`, `NGU`, and `NTU` tags (#reads supporting the four bases in normal) to the output. 3. Hybrid calling mode now outputs three VCFs, ClairS somatic variant calls, Clair3 normal germline variant calls, and Clair3 tumor germline variant calls. 4. Added the `--enable_clair3_germline_output` option to also output Clair3 normal germline variant calls, and Clair3 tumor germline variant calls (even when hybrid calling more is not enabled). Running time will increase by ~40%. 

*v0.1.0 (Mar 24)* : 1. Added support for Indel calling. ClairS Indel calling currently only supports ONT R10 data. To enable, use the `--enable_indel_calling` option. The Indel F1-score is ~73% with 50x/50x HCC1395/BL data. 2. Added an experimental `--normal_vcf_fn` to skip germline variant calling on normal BAM ([#7](https://github.com/HKU-BAL/ClairS/pull/7), contributor @[Xingyao](https://github.com/xingyaoc)). 3. Added `--hybrid_mode_vcf_fn` option to enable hybrid calling mode that combines de novo calling results and genotyping results without running the tool twice. Renamed the `--vcf_fn` to `--genotyping_mode_vcf_fn` for clarification. 4. Fixed a memory issue, memory consumption is now sub 100G for high coverage samples. 5. Fixed a conda environment issue in Singularity ([#3](https://github.com/HKU-BAL/ClairS/issues/3)). 6. Fixed zero division when no SNV was found ([#2](https://github.com/HKU-BAL/ClairS/issues/2), [#5](https://github.com/HKU-BAL/ClairS/issues/5)). 7. Added `AD` tag in the output.

*v0.0.1 (Jan 29, 2023)*: Initial release for early access.

---

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi [Revio](https://www.pacb.com/revio/) data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).
- Illumina NGS data as input, see [Illumina Quick Demo](docs/illumina_quick_demo.md).

### Quick start

After following [installation](#installation), you can run ClairS with one command:

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy
## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

------

## Pre-trained Models

ClairS trained both pileup and full-alignment models using GIAB samples, and carry on benchmarking on HCC1395-HCC1395BL pair dataset. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|        Platform         |       Model name       |    Chemistry /Instruments    |    Basecaller    |    Option (`-p/--platform`)    |   Reference   | Aligner  |
|:-----------------------:| :--------------------: | :--------------------------: | :----------------------: |:------------------------------:| :------: | ----------- |
|     ONT<sup>1</sup>     | r1041_e82_400bps_sup_v420 | R10.4.1, 5khz | Dorado SUP | `ont_r10_dorado_sup_5khz_ssrs` | GRCh38_no_alt | Minimap2 |
|     ONT<sup>1</sup>     | r1041_e82_400bps_sup_v420 | R10.4.1, 5khz | Dorado SUP |  `ont_r10_dorado_sup_5khz_ss`  | GRCh38_no_alt | Minimap2 |
|           ONT           | r1041_e82_400bps_sup_v410 | R10.4.1, 4khz | Dorado SUP |   `ont_r10_dorado_sup_4khz`    | GRCh38_no_alt | Minimap2 |
|           ONT           | r1041_e82_400bps_hac_v420 | R10.4.1, 5khz | Dorado HAC |   `ont_r10_dorado_hac_5khz`    | GRCh38_no_alt | Minimap2 |
|           ONT           | r1041_e82_400bps_hac_v410 | R10.4.1, 4khz | Dorado HAC |   `ont_r10_dorado_hac_4khz`    | GRCh38_no_alt | Minimap2 |
|           ONT           | r104_e81_sup_g5015 |        R10.4/R10.4.1, 4khz         |        Guppy5 SUP        |        `ont_r10_guppy`         | GRCh38_no_alt | Minimap2 |
|     ONT<sup>2</sup>     |  r941_prom_sup_g5014   |            R9.4.1, 4khz            |            Guppy5 SUP            |         `ont_r9_guppy`         | GRCh38_no_alt | Minimap2 |
|  Illumina<sup>1</sup>   |          ilmn          |        NovaSeq/HiseqX        |        -        |          `ilmn_ssrs`           |    GRCh38     | BWA-MEM  |
|  Illumina<sup>1</sup>   |          ilmn          |        NovaSeq/HiseqX        |        -        |           `ilmn_ss`            |    GRCh38     | BWA-MEM  |
| PacBio HiFi<sup>3</sup> |          hifi_sequel2          | Sequel II with Chemistry 2.0 | - |         `hifi_sequel2`         | GRCh38_no_alt | Minimap2 |
| PacBio HIFI<sup>1</sup> | hifi_revio | Revio with SMRTbell prep kit 3.0 | - |       `hifi_revio_ssrs`        | GRCh38_no_alt | Minimap2 |
| PacBio HIFI<sup>1</sup> | hifi_revio | Revio with SMRTbell prep kit 3.0 | - |        `hifi_revio_ss`         | GRCh38_no_alt | Minimap2 |

**Caveats <sup>1</sup>**: Starting from v0.4.0 version, ClairS will provide two model types. `ssrs` is a model trained initially with synthetic samples and then real samples augmented (e.g., `ont_r10_dorado_sup_5khz_ssrs`), `ss` is a model trained from synthetic samples (e.g., `ont_r10_dorado_sup_5khz_ss`). The `ssrs` model provides better performance and fits most usage scenarios. `ss` model can be used when missing a cancer-type in model training is a concern. In v0.4.0, four real cancer cell-line datasets (HCC1937, HCC1954, H1437, and H2009) covering two cancer types (breast cancer, lung cancer) published by [Park et al.](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1) were used for `ssrs` model training.

**Caveats <sup>2</sup>**: Although the r9(`r941_prom_sup_g5014`) model was trained on synthetic samples with r9.4.1 real data, the minimal AF cutoff, minimal coverage, and post-calling filtering parameters for the r9 model are copied from the r10 model, and are not optimized due to lack of real r9 data on a cancer sample with known truths.

**Caveats <sup>3</sup>**: The PacBio HiFi Sequel II model is experimental. It was trained but not tested with any real data with known truths. HG003 54x and HG004 52x were used, thus tumor depth coverage higher than 50x may suffer from lower recall rate. For testing, please downsample both tumor and normal to ~40x for the best performance of this experimental model.


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
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy, ont_r9_guppy, ilmn, hifi_sequel2, hifi_revio, etc}
  --output_dir ${OUTPUT_DIR}                   ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
mkdir -p ${OUTPUT_DIR}

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clairs:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clairs_latest.sif \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --normal_bam_fn ${INPUT_DIR}/normal.bam \    ## use your normal bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy, ont_r9_guppy, ilmn, hifi_sequel2, hifi_revio, etc}
  --output_dir ${OUTPUT_DIR} \                 ## output path prefix
  --conda_prefix /opt/conda/envs/clairs
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
mkdir -p ${CONDA_PREFIX}/bin/clairs_models
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
  --platform ${PLATFORM} \                   ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy, ont_r9_guppy, ilmn, hifi_sequel2, hifi_revio, etc}
  --output_dir ${OUTPUT_DIR}                 ## output path prefix
 
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
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options {ont_r10_dorado_sup_4khz, ont_r10_dorado_sup_5khz, ont_r10_dorado_hac_5khz, ont_r10_dorado_hac_4khz, ont_r10_guppy, ont_r9_guppy, ilmn, hifi_sequel2, hifi_revio}.
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
  -G GENOTYPING_MODE_VCF_FN, --genotyping_mode_vcf_fn GENOTYPING_MODE_VCF_FN
                        VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided.
  -H HYBRID_MODE_VCF_FN, --hybrid_mode_vcf_fn HYBRID_MODE_VCF_FN  
                        Enable hybrid calling mode that combines the de novo calling results and genotyping results at the positions in the VCF file given.
  --snv_min_qual SNV_MIN_QUAL                                                                                                                                                                              
                        If set, SNV variants with >SNV_MIN_QUAL will be marked as PASS, or LowQual otherwise.                                                             
  --indel_min_qual INDEL_MIN_QUAL                                                                                                                                                                          
                        If set, INDEL variants with >INDEL_MIN_QUAL will be marked as PASS, or LowQual otherwise.  
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
  --disable_phasing     EXPERIMENTAL: Disable phasing with longphase or whatshap. Usually leads to significant performance loss.
  --haplotagged_tumor_bam_provided_so_skip_intermediate_phasing_and_haplotagging
                        EXPERIMENTAL: Use haplotagged tumor bam as input, If enabled, will skip intermediate phasing and haplotagging, and the haplotype information
                        will be acquired from haplotagged tumor BAM. Default: disabled.
  --normal_vcf_fn NORMAL_VCF_FN
                        EXPERIMENTAL: Path to normal VCF file. Setting this will skip germline varaint calling on normal BAM file input.
  --enable_indel_calling
                        EXPERIMENTAL: Enable Indel calling, 'ont_r9_guppy' and 'ilmn' platforms are not supported. The calling time would increase significantly. default: disabled.
  --enable_clair3_germline_output
                        EXPERIMENTAL: Use Clair3 default calling settings than Clair3 fast calling setting for tumor and normal germline varaint calling. The calling time would increase ~40 percent, Default: disabled.
  --enable_verdict      EXPERIMENTAL: Use Verdict to tag the germline variant in CNV regions. We suggest using the parameter only for sample with tumor purity lower than 0.8, Default: disabled
  --use_heterozygous_snp_in_normal_sample_and_tumor_bam_for_intermediate_phasing USE_HETEROZYGOUS_SNP_IN_NORMAL_SAMPLE_AND_TUMOR_BAM_FOR_INTERMEDIATE_PHASING
                        EXPERIMENTAL: Use the heterozygous SNPs in normal VCF called by Clair3 and the tumor BAM for intermediate phasing. Option: {True, False}.
                        Default: True.
  --use_heterozygous_snp_in_tumor_sample_and_tumor_bam_for_intermediate_phasing USE_HETEROZYGOUS_SNP_IN_TUMOR_SAMPLE_AND_TUMOR_BAM_FOR_INTERMEDIATE_PHASING
                        EXPERIMENTAL: Use the heterozygous SNPs in tumor VCF called by Clair3 and the tumor BAM for intermediate phasing. Option: {True, False}.
                        Default: False.
  --use_heterozygous_snp_in_normal_sample_and_normal_bam_for_intermediate_phasing USE_HETEROZYGOUS_SNP_IN_NORMAL_SAMPLE_AND_NORMAL_BAM_FOR_INTERMEDIATE_PHASING
                        EXPERIMENTAL: Use the heterozygous SNPs in normal VCF called by Clair3 and the normal BAM for intermediate phasing. Option: {True, False}.
                        Default: False.
  --use_heterozygous_snp_in_tumor_sample_and_normal_bam_for_intermediate_phasing USE_HETEROZYGOUS_SNP_IN_TUMOR_SAMPLE_AND_NORMAL_BAM_FOR_INTERMEDIATE_PHASING
                        EXPERIMENTAL: Use the heterozygous SNPs in tumor VCF called by Clair3 and the normal BAM for intermediate phasing. Option: {True, False}.
                        Default: False.
  --use_heterozygous_indel_for_intermediate_phasing USE_HETEROZYGOUS_INDEL_FOR_INTERMEDIATE_PHASING
                        EXPERIMENTAL: Use the heterozygous Indels in normal and tumor VCFs called by Clair3 for intermediate phasing. Option: {True, False}. Default: True.
  --use_longphase_for_intermediate_haplotagging USE_LONGPHASE_FOR_INTERMEDIATE_HAPLOTAGGING
                        EXPERIMENTAL: Use the longphase instead of whatshap for intermediate haplotagging. Option: {True, False}. Default: True.
  --indel_output_prefix INDEL_OUTPUT_PREFIX
                        Prefix for Indel output VCF filename. Default: indel.
  --indel_pileup_model_path INDEL_PILEUP_MODEL_PATH
                        Specify the path to your own somatic calling indel pileup model.
  --indel_full_alignment_model_path INDEL_FULL_ALIGNMENT_MODEL_PATH
                        Specify the path to your own somatic calling indel full-alignment model.
```

#### Call SNVs in one or mutiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -C chr21,chr22
```

#### Call SNVs in one specific region using the `-r/--region` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -r chr20:1000000-2000000
```

#### Call SNVs at interested variant sites (genotyping) using the `-G/--genotyping_mode_vcf_fn` parameter

```bash
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -G input.vcf
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
./run_clairs -T tumor.bam -N normal.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -B input.bed
```

------


## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
