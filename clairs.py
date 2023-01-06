import sys
from importlib import import_module

DATA_PREP_SCRIPTS_FOLDER="src"
DEEP_LEARNING_FOLDER="clairs"

deep_learning_folder = [
    "train",
    "predict",
    "call_variants",

    "CallVarBam",
    "CallVariants",
    "train",
    "Predict",
    "Predict_with_logits",
    "Train_multi-task",
    "Predict_multi_class",
    "Train_multi-task2",
    "Predict_multi_class_add",
    "Train_torch",
    "predict",
    "call_variants",
    "call_variants_from_bam",
    "Predict_all",
    'train_adjacent',
    'predict_adjacent',
    'call_variants_adjacent',
    "predict_cal_entropy",
    'train_tf',
    "predict_tf"
]

REPO_NAME = "clair-somatic"
data_preprocess_folder = [
    "get_candidates",
    "create_tensor",
    "create_bin",
    "split_bam",
    "mix_chunk_bam",
    "extract_candidates",
    "filter_reference_calls",
    "create_tensor_pileup",
    'find_tumor_truth_in_normal',
    "update_variant",
    "create_pair_tensor",
    "sort_vcf",
    "compare_vcf",
    "create_pair_tensor_pileup",
    "extract_pair_candidates",
    'realign_reads',
    'realign_variants',
    'select_hetero_snp_for_phasing',
    "gen_contaminated_bam",
    "haplotype_filtering",
    "merge_vcf",
    'cal_af_distribution',
    'add_back_missing_variants_in_genotyping',

    "get_candidates",
    "extract_candidates_test",
    "create_tensor",
    "create_bin",
    "split_bam",
    "PlotAlignment",
    "mix_chunk_bam",
    "extract_candidates",
    "Vcf2Bed",
    "PlotAF",
    "filter_reference_calls",
    "create_tensor_adjacent",
    "create_tensor_pileup",
    'find_tumor_truth_in_normal',
    "update_variant",
    "create_pair_tensor",
    "sort_vcf",
    "compare_vcf",
    "extract_candidates_with_af",
    "extract_candidates_with_af_ru",
    'create_tensor_test',
    "create_pair_tensor_adjacent",
    "create_pair_tensor_pileup",
    "create_bin_adjacent",
    'PlotAlignment_no_chr',
    "extract_pair_candidates",
    'create_tensor_output_read_name_for_seqc_pb',
    "compare_vcf_appy_qual_filter",
    'realign_reads',
    'realign_variants',
    "realign_variants_bak",
    "compare_vcf_appy_qual_strva_v2",
    'select_hetero_snp_for_phasing',
    "create_pair_tensor_no_low_bq_read",
    "gen_contaminated_bam",
    "haplotype_filtering",
    "merge_vcf",
    "create_pair_tensor_pileup_rm_rn",
    "extract_pair_candidates_rm_rn"
]


def directory_for(submodule_name):
    if submodule_name in deep_learning_folder:
        return DEEP_LEARNING_FOLDER
    if submodule_name in data_preprocess_folder:
        return DATA_PREP_SCRIPTS_FOLDER
    return ""


def print_help_messages():
    from textwrap import dedent
    print(dedent("""\
        {0} submodule invocator:
            Usage: python clair-somatic.py [submodule] [options of the submodule]
        Available data preparation submodules:\n{1}
        Available clair submodules:\n{2}
        """.format(
            REPO_NAME,
            "\n".join("          - %s" % submodule_name for submodule_name in data_preprocess_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in deep_learning_folder),
        )
    ))


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print_help_messages()
        sys.exit(0)

    submodule_name = sys.argv[1]
    if (
        submodule_name not in deep_learning_folder and
        submodule_name not in data_preprocess_folder
    ):
        sys.exit("[ERROR] Submodule %s not found." % (submodule_name))

    directory = directory_for(submodule_name)
    submodule = import_module("%s.%s" % (directory, submodule_name))

    sys.argv = sys.argv[1:]
    sys.argv[0] += (".py")

    # Note: need to make sure every submodule contains main() method
    submodule.main()

    sys.exit(0)


if __name__ == "__main__":
    main()