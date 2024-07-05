from argparse import ArgumentParser

import numpy as np
from scipy.interpolate import BSpline
from sklearn.linear_model import LinearRegression


def create_bspline_basis(x, df, degree=3):
    """Create B-spline basis for given data"""
    n_knots = df - degree + 1
    knots = np.linspace(np.min(x), np.max(x), n_knots)
    knots = np.concatenate(([knots[0]] * degree, knots, [knots[-1]] * degree))
    spline = BSpline(knots, np.eye(len(knots) - degree - 1), degree)
    return np.row_stack([spline(xi) for xi in x])


def correctLogR(tumor_logr_file, gc_content_file, replication_timing_file, tumor_logr_correction_output_file, sample_name):
    tumor_logr = open(tumor_logr_file, 'r')
    gc_content = open(gc_content_file, 'r')
    replication_timing = open(replication_timing_file, 'r')
    tumor_logr_dict = dict()
    gc_content_dict = dict()
    replication_timing_dict = dict()
    for idx, tumor_logr_line in enumerate(tumor_logr.readlines()):
        if idx == 0:
            continue
        tumor_logr_info = tumor_logr_line.strip().split('\t')
        chr = tumor_logr_info[0]
        pos = tumor_logr_info[1]
        logr = tumor_logr_info[2]
        key = (str(chr), str(pos))
        tumor_logr_dict[key] = logr
    for idx, gc_content_line in enumerate(gc_content.readlines()):
        if idx == 0:
            continue
        gc_content_info = gc_content_line.strip().split('\t')
        chr = 'chr' + str(gc_content_info[1])
        pos = gc_content_info[2]
        gccontent = ('\t').join(gc_content_info[3:])
        key = (str(chr), str(pos))
        gc_content_dict[key] = gccontent
    for idx, replication_timing_line in enumerate(replication_timing.readlines()):
        if idx == 0:
            continue
        replication_timing_info = replication_timing_line.strip().split('\t')
        chr = 'chr' + str(replication_timing_info[1])
        pos = replication_timing_info[2]
        replicationtiming = ('\t').join(replication_timing_info[3:])
        key = (str(chr), str(pos))
        replication_timing_dict[key] = replicationtiming

    overlap_keys = tumor_logr_dict.keys()
    tumor_logr_dict_overlap = {key: tumor_logr_dict[key] for key in overlap_keys}
    gc_content_dict_overlap = {key: gc_content_dict[key] for key in overlap_keys}

    tumor_logr_dict_values = np.array(list(tumor_logr_dict_overlap.values())).astype(float)
    gc_content_dict_values = [gc_content_value.split('\t') for gc_content_value in list(gc_content_dict_overlap.values())]
    gc_content_dict_values = np.array(gc_content_dict_values).astype(float)

    corr_gc = np.abs(np.corrcoef(gc_content_dict_values, tumor_logr_dict_values, rowvar=False))[0, 1:]

    index_1kb = 5
    index_max = 11
    maxGCcol_insert = np.argmax(corr_gc[:(index_1kb + 1)])
    maxGCcol_amplic = np.argmax(corr_gc[(index_1kb + 2):(index_max + 1)]) + (index_1kb + 2)

    replication_timing_dict_overlap = {key: replication_timing_dict[key] for key in overlap_keys}
    replication_timing_dict_values = [replication_timing_value.split('\t') for replication_timing_value in
                              list(replication_timing_dict_overlap.values())]
    replication_timing_dict_values = np.array(replication_timing_dict_values).astype(float)

    corr_rep = np.abs(np.corrcoef(replication_timing_dict_values, tumor_logr_dict_values, rowvar=False))[0, 1:]

    maxreplic = np.argmax(corr_rep)

    GC_insert_bspline = create_bspline_basis(gc_content_dict_values[:, maxGCcol_insert], df=5)
    GC_amplic_bspline = create_bspline_basis(gc_content_dict_values[:, maxGCcol_amplic], df=5)
    replic_bspline = create_bspline_basis(replication_timing_dict_values[:, maxreplic], df=5)

    X = np.hstack([GC_insert_bspline, GC_amplic_bspline, replic_bspline])
    y = tumor_logr_dict_values.reshape(-1, 1)
    model = LinearRegression(fit_intercept=True).fit(X, y)
    residuals = y - model.predict(X)
    tumor_logr_dict_values_after = residuals.flatten()

    tumor_logr_dict_overlap_after = {key: tumor_logr_dict_values_after[idx] for idx, key in enumerate(overlap_keys)}

    output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
    tumor_logr_correction_output = open(tumor_logr_correction_output_file, 'w')
    tumor_logr_correction_output.write(output_header)
    for key, value in tumor_logr_dict_overlap_after.items():
        tumor_logr_correction_string = '\t'.join(map(str, key)) + '\t' + str(value) + '\n'
        tumor_logr_correction_output.write(tumor_logr_correction_string)


def main():
    parser = ArgumentParser(description="Correct Tumor Sample LogR")

    parser.add_argument('--tumor_logr_file', type=str,
                        default=None,
                        help="Path of tumor sample LogR")

    parser.add_argument('--gc_content_file', type=str,
                        default=None,
                        help="Path of 1kG GC content file")

    parser.add_argument('--replication_timing_file', type=str,
                        default=None,
                        help="Path of 1kG replication timing file")

    parser.add_argument('--tumor_logr_correction_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample corrected LogR")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")


    global args
    args = parser.parse_args()

    correctLogR(args.tumor_logr_file, args.gc_content_file, args.replication_timing_file, args.tumor_logr_correction_output_file, args.sample_name)


if __name__ == "__main__":
    main()