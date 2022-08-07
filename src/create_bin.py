import sys
import logging
from argparse import ArgumentParser, SUPPRESS
from shared.utils import str2bool

logging.basicConfig(format='%(message)s', level=logging.INFO)

def Run(args):

    if args.normal_tensor_fn is None:
        import clair_somatic.utils_pair as utils
    elif args.pileup:
        import clair_somatic.utils_pileup as utils
    else:
        import clair_somatic.utils as utils
    utils.setup_environment()
    logging.info("Loading the dataset ...")
    
    utils.get_training_array(
        normal_tensor_fn=args.normal_tensor_fn,
        tumor_tensor_fn=args.tumor_tensor_fn,
        var_fn=args.var_fn,
        bed_fn=args.bed_fn,
        bin_fn=args.bin_fn,
        shuffle=args.shuffle,
        is_allow_duplicate_chr_pos=args.allow_duplicate_chr_pos,
        chunk_id=args.chunk_id-1 if args.chunk_id else None, # 1-base to 0-base
        chunk_num=args.chunk_num,
        pileup=args.pileup,
        platform=args.platform,
        merge_bins=args.merge_bins)
    logging.info("Finish!")


def main():
    parser = ArgumentParser(description="Combine the variant and non-variant tensors and convert them to a binary")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--normal_tensor_fn', type=str, default=None,
                        help="Normal tensor input, required")

    parser.add_argument('--tumor_tensor_fn', type=str, default=None,
                        help="Tumor tensor input, required")

    parser.add_argument('--var_fn', type=str, default=None,
                        help="Truth variants list input, required")

    parser.add_argument('--bin_fn', type=str, default=None,
                        help="Output a binary tensor file, required")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="High confident genome regions input in the BED format")

    parser.add_argument('--shuffle', action='store_true',
                        help="Shuffle the inputs")

    parser.add_argument('--allow_duplicate_chr_pos', action='store_true',
                        help="Allow duplicated chromosome:position in the tensor input")

    # options for internal process control
    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Maximum non-variant ratio against variant in the training data
    parser.add_argument('--maximum_non_variant_ratio', type=float, default=None,
                        help=SUPPRESS)

    ## Path to the variant candidate details
    parser.add_argument('--candidate_details_fn_prefix', type=str, default=None,
                        help=SUPPRESS)

    ## Path to the variant candidate details
    parser.add_argument('--merge_bins', type=str2bool, default=False,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
