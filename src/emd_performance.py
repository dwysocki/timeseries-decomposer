#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import defaultdict
from itertools import chain, tee
from os import makedirs, path
from os.path import isdir
from subprocess import check_output
from sys import stdout

import numpy as np
import matplotlib.pyplot as plt

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
        default=10,
        help="")
    parser.add_argument("--samples-per-bin", type=int,
        default=10,
        help="")
    parser.add_argument("--seed", type=int,
        help="Random number generator seed.")

    parser.add_argument("--num-cycles", type=float,
        default=2.0,
        help="")
    parser.add_argument("--time-step", type=float,
        default=0.01,
        help="")
    parser.add_argument("--r-seed", type=int,
        default=0,
        help="Random number generator seed for R script.")

    args = parser.parse_args()

    if args.output is not None:
        make_sure_path_exists(args.output)

    return args


def main():
    args = get_args()

    # separators values for standard deviation bins
    bins = np.linspace(0.0, args.noise_sd_max, args.num_bins+1,
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
    #  statM: [mean(statM from bin1), ..., mean(statM from binN)]}
    stats = compressed_average_dicts(stat_samples)

    # print stats table to stdout
    output_table = np.column_stack(chain([bins[:-1], bins[1:]], stats.values()))
    header = "\t".join(chain(["lo", "hi"], stats.keys()))
    np.savetxt(stdout.buffer, output_table, header=header, delimiter="\t")

    if args.output is not None:
        make_plots(bins, stats, args)



def run_emd(args, sd):
    # create calling sequence
    call = ["./emd_concept.r",
            # positional arguments
            args.signal_type, args.noise_type, args.criteria,
            # keyword arguments
            "--num-cycles", args.num_cycles,
            "--time-step", args.time_step,
            "--seed", args.r_seed,
            # most importantly, give the standard deviation for this call
            "--noise-sd", sd]
    # ensure all arguments are of type str
    call = map(str, call)
    results = check_output(call)

    return parse_emd(results)


def parse_emd(results):
    header_str, stats_str = results.splitlines()

    header = map(lambda b: b.decode(), header_str.strip().split())
    stats  = map(float, stats_str.strip().split())

    return dict(zip(header, stats))


def make_plots(bins, stats, args):
    for stat_name, stat_vals in stats.items():
        plot_histogram(bins, stat_vals, stat_name, args.output)


def plot_histogram(bins, values, name, output):
    left = bins[:-1]
    right = bins[1:]
    width = right - left

    fig, ax = plt.subplots()

    ax.bar(left=left, width=width, height=values)

    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel(name)

    fig.savefig(path.join(output, name+"_histogram"))


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


def make_sure_path_exists(path):
    """make_sure_path_exists(path)

    Creates the supplied *path* if it does not exist.
    Raises *OSError* if the *path* cannot be created.

    **Parameters**

    path : str
        Path to create.

    **Returns**

    None
    """
    try:
        makedirs(path)
    except OSError:
        if not isdir(path):
            raise


if __name__ == "__main__":
    exit(main())
