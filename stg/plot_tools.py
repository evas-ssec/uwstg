#!/usr/bin/env python
# encoding: utf-8
"""
This module draws simple plots for space and time gridded data in flat binary
files.

:author:       Eva Schiffer (evas)
:contact:      eva.schiffer@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
:date:         Jan 2014
:license:      GNU GPLv3

Copyright (C) 2014 Space Science and Engineering Center (SSEC),
 University of Wisconsin-Madison.
"""
__docformat__ = "restructuredtext en"

import os
from glob import glob
import numpy

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import keoni.fbf.workspace as Workspace

DEFAULT_FILE_PATTERN = "*.real4.*.*"
DEFAULT_FILL_VALUE   = numpy.nan
DEFAULT_DPI          = 150

def plot_binary(bf, in_dir='.',
                fill_value=DEFAULT_FILL_VALUE,
                dpi_to_use=DEFAULT_DPI,
                vmin=None, vmax=None):
    
    # get the data from the file
    var_workspace = Workspace.Workspace(dir=in_dir)
    fbf_attr_name = bf.split(".")[0]
    raw_data      = var_workspace[fbf_attr_name][:]
    # TODO, there's probably a better way to show 3D data
    raw_data   = numpy.nansum(raw_data, axis=0) if len(raw_data.shape) > 2 else raw_data
    
    
    # mask the data based on the fill value
    masked_data   = numpy.ma.masked_where(raw_data == fill_value, raw_data)
    print masked_data.min(), masked_data.max()
    print masked_data.shape
    
    # plot the figure
    plt.figure()
    plt.imshow(masked_data, vmin=vmin, vmax=vmax)
    plt.bone()
    plt.colorbar()
    plt.savefig("plot_binary.%s.png" % fbf_attr_name, dpi=dpi_to_use)
    plt.close()

def sci_float(x):
    x = x.replace("\"", "")
    x = x.replace("\'", "")
    return float(str(x))

def main():
    from argparse import ArgumentParser
    description = """
Plot binary files using matplotlib.
    """
    parser = ArgumentParser(description=description)
    parser.add_argument("-f", dest="fill_value", default=DEFAULT_FILL_VALUE, type=sci_float,
            help="Specify the fill_value of the input file(s)")
    parser.add_argument('--vmin', dest="vmin", default=None, type=int,
            help="Specify minimum brightness value. Defaults to minimum value of data.")
    parser.add_argument('--vmax', dest="vmax", default=None, type=int,
            help="Specify maximum brightness value. Defaults to maximum value of data.")
    parser.add_argument("-p", dest="pattern",
            help="filename pattern to search the current directory for")
    parser.add_argument("binary_files", nargs="*",
            help="list of flat binary files to be plotted in the current directory")
    parser.add_argument('-d', '--dpi', dest="dpi", default=DEFAULT_DPI, type=float,
            help="Specify the dpi for the resulting figure, higher dpi will result in larger figures and longer run times")
    args = parser.parse_args()

    workspace = '.'
    binary_files = args.binary_files
    if not args.binary_files and not args.pattern:
        args.pattern = DEFAULT_FILE_PATTERN
    if args.pattern:
        workspace = os.path.split(args.pattern)[0]
        binary_files = [ os.path.split(x)[1] for x in glob(args.pattern) ]

    for bf in binary_files:
        print "Plotting '%s'" % (bf,)
        try:
            plot_binary(bf, in_dir='.',
                        fill_value=args.fill_value,
                        dpi_to_use=args.dpi,
                        vmin=args.vmin, vmax=args.vmax)
        except StandardError as e:
            print "Could not plot '%s'" % (bf,)
            if hasattr(e, "msg"): print e,e.msg
            else: print e

if __name__ == "__main__":
    import sys

    sys.exit(main())

