#!/usr/bin/env python

import argparse
import plotext
import numpy as np
import re
import plotext as plt

def parse_cols(file, columns=[]):
    with open(file, 'r') as f:
        lines = list(filter(lambda x: x.startswith('@ s'), f.read().splitlines()))
    indices = list(map(lambda x: int(x.lstrip('@ ').split()[0].lstrip('s')), lines))
    labels = [re.findall(r"\"([A-Za-z0-9_]+)\"", s)[0] for s in lines]
    if not columns:
        return {k: v for k, v in zip(indices, labels)}
    else:
        return {k: v for k, v in zip(indices, labels) if v in list(map(lambda x: x.lower(), columns))}


def get_axis_labels(file):
    with open(file, 'r') as f:
        lines = list(filter(lambda x: 'axis  label' in x, f.read().splitlines()))
    labels = [s.split('"')[-2] for s in lines]
    assert len(labels) == 2
    return labels


def get_title(file):
    with open(file, 'r') as f:
        lines = list(filter(lambda x: 'title' in x, f.read().splitlines()))
    titles = [s.split('"')[-2] for s in lines]
    assert len(titles) == 1
    return titles[0]


def main(file, columns=[]):
    usecols = parse_cols(file, columns=columns)
    data = np.loadtxt(file, comments=['#', '@'], usecols=[0] + [i + 1 for i in usecols.keys()])
    labels = get_axis_labels(file)
    title = get_title(file)

    for i, label in usecols.items():
        plt.plot(data[:,i + 1], label=label)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="plot_xvg",
        usage="plot_xvg.py [file.xvg] -cols",
        description="Plots a gromacs .xvg file to the terminal using plotext."
    )
    parser.add_argument(
        dest="xvg_file",
        help="The xvg file containing the data."
    )
    parser.add_argument(
        "-cols",
        "--columns",
        required=False,
        default='all',
        dest="columns",
        help="""The fields that should be plotted. Can either be a list of [ener\
                gy, potential], or all."""
    )
    args = parser.parse_args()
    if args.columns == 'all':
        args.columns = {}
    main(args.xvg_file, args.columns)

