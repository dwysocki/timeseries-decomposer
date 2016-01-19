#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import defaultdict
from itertools import tee
from os import path
import numpy as np

def get_args():
    parser = ArgumentParser()

    parser.add_argument("signal_type",
        help="")
    parser.add_argument("noise_type",
        help="")
    parser.add_argument("criteria", type=int,
        help="")

    parser.add_argument("-o", "--output",
        help="Output directory for plots, if desired.")

    parser.add_argument("--noise-sd-max", type=float,
        default=1.0,
        help="")
    parser.add_argument("--num-bins", type=int,
        default=11)
    parser.add_argument("--samples-per-bin", type=int,
        default=10)

    args = parser.parse_args()

    return args


def main():
    args = get_args()

    # separators values for standard deviation bins
    bins = np.linspace(0.0, args.noise_sd_max, args.num_bins,
                       endpoint=True)
    # a list of binned samples of noise standard deviations to use
    sd_samples = [np.random.uniform(a, b, args.samples_per_bin)
                  for a, b
                  in pairwise(bins)]
    # a list of computed statistics for each binned statistic
    stat_samples = [[run_emd(args, sd) for sd in sd_bin]
                    for sd_bin in sd_samples]
    # average each bin and create a map
    # {stat1: [mean(stat1 from bin1), ..., mean(stat1 from binN)],
    #  ...,
    #  statM; [mean(statM from bin1), ..., mean(statM from binN)]}
    stats = compressed_average_dicts(stat_samples)

    for stat_name, stat_vals in stats.items():
        print("#### {stat_name} ####".format(stat_name=stat_name))

        for (lo, hi), stat in zip(pairwise(bins), stat_vals):
            print("{lo:.3f} <= x < {hi:.3f} = {stat:.3f}"
                  .format(lo=lo, hi=hi, stat=stat))

        print()


def run_emd(args, sd):
    return {"MSE": np.random.uniform(0, 1),
            "R^2": np.random.uniform(10, 11)}


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def compress_dicts(list_of_dicts):
    """
    Takes a list of dictionaries
    [{k_1: v_11, ..., k_N: v_N1}, ..., {k_1: v_1M, ..., k_N: v_NM}]
    and compresses it into a single dictionary of lists
    {k_1: [v_11, ..., v_1M], ..., k_N: [v_N1, ..., v_NM]}
    """
    # initialize dict
    ret = defaultdict(list)

    for d in list_of_dicts:
        for k, v in d.items():
            ret[k].append(v)

    return ret


def average_dict(dict_of_lists, copy=True):
    """
    Takes a dictionary whose values are numeric lists, and returns a dictionary
    with the same keys, and the values averaged. If copy=False, updates the
    dictionary in-place.
    """
    # initialize dict if copy=True, else use same dict
    ret = {k: None for k in dict_of_lists} if copy else dict_of_lists

    for k, v in dict_of_lists.items():
        ret[k] = np.mean(v)

    return ret


def compressed_average_dicts(list_of_list_of_dicts):
    """
    Takes a 2D list of dictionaries, and returns a dictionary mapping
    keys from the original nested dicts to a list of the average from each row.

    In other words, takes a list of list of dictionaries:
    [[{k_1: v_111, ..., k_N: v_i11}, ..., {k_1: v_1j1, ..., k_N: v_ij1}],
     ...,
     [{k_1: v_11k, ..., k_N: v_i1k}, ..., {k_1: v_1jk, ..., k_N: v_ijk}]]
    and returns a dictionary of lists
    {k_1: [mean(v_111, ..., v_1j1), ..., mean(v_11k, ..., v_1jk)],
     ...,
     k_N: [mean(v_i11, ..., v_ij1), ..., mean(v_i1k, ..., v_ijk)]}
    """
    list_of_dicts_of_lists = [compress_dicts(list_of_dicts)
                              for list_of_dicts in list_of_list_of_dicts]
    list_of_dicts_of_averages = [average_dict(x, copy=False)
                                 for x in list_of_dicts_of_lists]
    dict_of_averages = compress_dicts(list_of_dicts_of_averages)

    return dict_of_averages


if __name__ == "__main__":
    exit(main())
