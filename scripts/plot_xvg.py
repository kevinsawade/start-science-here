#!/usr/bin/env python

import argparse
import plotext
import numpy as np
import re
import plotext as plt

def parse_cols(file, columns=[]):
    with open(file, 'r') as f:
        lines = f.read()
    energy = 'gmx energy' in lines
    lines = list(filter(lambda x: x.startswith('@ s'), lines.splitlines()))
    if energy:
        indices = list(map(lambda x: int(x.lstrip('@ ').split()[0].lstrip('s')), lines))
    else:
        indices = [0]
    labels = [s.split('"')[-2] for s in lines]
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
    if len(titles) == 0:
        raise Exception(f"Could not find title in {file}")
    elif len(titles) == 1:
        return titles[0]
    else:
        return ' '.join(titles)


def main(files, columns=[], running_average=0, grid=False, plot_labels=[]):
    lbl = 0
    for file in files:
        usecols = parse_cols(file, columns=columns)
        data = np.loadtxt(file, comments=['#', '@'], usecols=[0] + [i + 1 for i in usecols.keys()])
        labels = get_axis_labels(file)
        title = get_title(file)

        assert usecols

        for i, label in usecols.items():
            if labels:
                try:
                    label = plot_labels[lbl]
                except IndexError:
                    pass
            plt.plot(data[:,i + 1], label=label)
            if running_average:
                y = data[:, i + 1]
                y = np.convolve(y, np.ones(running_average)/running_average, mode='valid')
                x = np.arange(len(y)) + running_average/2
                plt.plot(x, y, label=label + f' running avg over {running_average}', color='red')
            lbl += 1
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    if grid:
        plt.grid(True, True)
    plt.show()

    _ = input("Press any button to clear plot")
    plt.clt()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="plot_xvg",
        usage="plot_xvg.py [file.xvg] -cols",
        description="Plots a gromacs .xvg file to the terminal using plotext."
    )
    parser.add_argument(
        dest="xvg_file",
        nargs='+',
        help="The xvg file containing the data."
    )
    parser.add_argument(
        "-cols",
        "--columns",
        required=False,
        default='all',
        dest="columns",
        help="""The fields that should be plotted. Can either be a list of [energy, potential], or all."""
    )
    parser.add_argument(
        "-rav",
        "--running-avg",
        required=False,
        default=0,
        type=int,
        dest='running_average',
        help="""Adds n steps running averages to the plots. For example: -rav 10 sets a 10 step runnign average for all plots."""
    )
    parser.add_argument(
        "-grid",
        required=False,
        default=False,
        action='store_true',
        dest='grid',
        help="Whether to add a grid to the figure. This option does not take arguments. Set to set grid true."
    )
    parser.add_argument(
        '-lbl',
        '--labels',
        nargs='+',
        required=False,
        default = [],
        dest = 'labels',
        help = 'Manually set plot labels. This option takes as many labels as there are plots to be plotted. Further provided labels will be discarded.'
    )
    args = parser.parse_args()
    if args.columns == 'all':
        args.columns = {}
    main(args.xvg_file, args.columns, args.running_average, args.grid, args.labels)
