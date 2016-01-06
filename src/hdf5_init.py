"""
Creates an HDF5 file whose groups mimic the subdirectory structure of a given
directory.
"""

from os import walk
from os.path import relpath
from sys import argv, stderr
import h5py

def empty_directories(top):
    return (
        relpath(dirpath, top)
        for dirpath, dirs, files
        in walk(top)
        if not dirs
    )


def main(input_dir, hdf5_fname):
    with h5py.File(hdf5_fname, "w") as hdf5_file:
        for d in empty_directories(input_dir):
            print("Creating group {}".format(d),
                  file=stderr)
            hdf5_file.create_group(d)


if __name__ == "__main__":
    exit(main(*argv[1:]))
