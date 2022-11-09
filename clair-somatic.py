import sys
from importlib import import_module

DATA_PREP_SCRIPTS_FOLDER="src"
DEEP_LEARNING_FOLDER="clair_somatic"
POST_PROCESS_SCRIPTS_FOLDER="clair_somatic.metrics"

deep_learning_folder = [
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
    "Predict_all"
]

REPO_NAME = "clair-somatic"
data_preprocess_folder = [
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
    "extract_candidates_with_af_ru"
]

post_process_scripts_folder = [
    'GetOverallMetrics',
]

def directory_for(submodule_name):
    if submodule_name in deep_learning_folder:
        return DEEP_LEARNING_FOLDER
    if submodule_name in data_preprocess_folder:
        return DATA_PREP_SCRIPTS_FOLDER
    if submodule_name in post_process_scripts_folder:
        return POST_PROCESS_SCRIPTS_FOLDER
    return ""


def print_help_messages():
    from textwrap import dedent
    print(dedent("""\
        {0} submodule invocator:
            Usage: python clair3.py [submodule] [options of the submodule]
        Available data preparation submodules:\n{1}
        Available clair submodules:\n{2}
        Available post processing submodules:\n{3}
        """.format(
            REPO_NAME,
            "\n".join("          - %s" % submodule_name for submodule_name in data_preprocess_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in deep_learning_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in post_process_scripts_folder),
        )
    ))


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print_help_messages()
        sys.exit(0)

    submodule_name = sys.argv[1]
    if (
        submodule_name not in deep_learning_folder and
        submodule_name not in data_preprocess_folder and
        submodule_name not in post_process_scripts_folder
    ):
        sys.exit("[ERROR] Submodule %s not found." % (submodule_name))

    directory = directory_for(submodule_name)
    submodule = import_module("%s.%s" % (directory, submodule_name))

    # filter arguments (i.e. filter clair3.py) and add ".py" for that submodule
    sys.argv = sys.argv[1:]
    sys.argv[0] += (".py")

    # Note: need to make sure every submodule contains main() method
    submodule.main()

    sys.exit(0)


if __name__ == "__main__":
    main()
