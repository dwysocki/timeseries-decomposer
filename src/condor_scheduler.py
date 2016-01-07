from argparse import ArgumentParser
import numpy as np
import itertools as it
from tempfile import mkdtemp, NamedTemporaryFile

from os import walk
from sys import stdin, stdout, stderr



_scheduler_prefix = "condor_scheduler-"


def get_args():
    parser = ArgumentParser()

    parser.add_argument("-c", "--command", type=str,
        help="Command to run on each line of stdin. Follows xargs syntax.")
    parser.add_argument("-s", "--min-size", type=int,
        help="Minimum number of lines to process in a single condor job.")
    parser.add_argument("-g", "--max-jobs", type=int, default=None,
        help="Maximum number of condor jobs to schedule (optional).")
    parser.add_argument("-e", "--executable", type=str, default="/bin/bash",
        help="Executable to run command with (default /bin/bash).")
    parser.add_argument("--log", type=str, default="condor.log",
        help="")
    parser.add_argument("--err", type=str, default="condor.err",
        help="")
    ## files which are created as temp files if not specified
    ## TODO: explain this in the help text

    parser.add_argument("--condor-file", type=str, default=None,
        help="")
    parser.add_argument("--dag-file", type=str, default=None,
        help="")
    parser.add_argument("--script-file", type=str, default=None,
        help="")
    parser.add_argument("--input-dir", type=str, default=None,
        help="")

    args = parser.parse_args()

    return args


def main():
    args = get_args()

    # read stdin into a list of its lines
    lines = stdin.readlines()
    # partition the lines into groups
    groups = grouper(lines, args.min_size, args.max_jobs)
    #
    script = make_script(args.script_file, args.command, args.executable)
    #
    inputs = make_input_files(args.input_dir, groups)
    #
    condor = make_condor(args.condor_file, script, args.log, args.err)
    #
    dag    = make_dag(args.dag_file, condor, inputs)


def open_tmpfile(fname, suffix):
    if fname is None:
        return NamedTemporaryFile("w",
                                  prefix=_scheduler_prefix, suffix=suffix,
                                  delete=False)
    else:
        return open(fname, "w")


def write_to_tmpfile(fname, suffix, contents):
    with open_tmpfile(fname, suffix) as f:
        f.write(contents)
        return f.name


script_template = """\
#!{}

xargs {}
"""

def make_script(fname, command, executable):
    return write_to_tmpfile(fname, "",
                            script_template.format(executable, command))


def make_input_file(dir, group):
    with NamedTemporaryFile(mode="w",
                            prefix=_scheduler_prefix, suffix=".input",
                            dir=dir, delete=False) as f:
        # write group to file
        f.writelines([line for line in group if line is not None])

        # return path to file
        return f.name


def make_input_files(dir, groups):
    if dir is None:
        dir = mkdtemp(prefix=_scheduler_prefix)

    return [make_input_file(dir, group) for group in groups]


condor_template = """\
executable = {}
input = $(input)
log = {}
error = {}
queue
"""


def make_condor(fname, script, log, err):
    return write_to_tmpfile(fname, ".condor",
                            condor_template.format(script, log, err))


dag_template = "JOB {0} {1}\nVARS {0} input={0}\n"


def make_dag(fname, condor, inputs):
    with open_tmpfile(fname, ".dag") as f:
        for input in inputs:
            f.write(dag_template.format(input, condor))


def partition(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks

    >>> partition('ABCDEFG', 3, 'x')
    'ABC' 'DEF' 'Gxx'
    """
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)


def grouper(iterable, min_size, max_groups=None):
    group_size = min_size \
        if max_groups is None \
        else int(np.ceil(len(iterable) / max_groups))

    return partition(iterable, group_size)


if __name__ == "__main__":
    exit(main())
