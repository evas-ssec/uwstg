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
import matplotlib.cm          as cm

from stg.stg_util  import make_lat_lon_grids
from stg.constants import *

import keoni.fbf.workspace as Workspace
from mpl_toolkits.basemap import Basemap

DEFAULT_FILE_PATTERN = "*.real4.*.*"
DEFAULT_FILL_VALUE   = numpy.nan
DEFAULT_DPI          = 150
DEFAULT_LEVELS_NUM   = 50
DEFAULT_LON_RANGE    = [-180, 180]
DEFAULT_LAT_RANGE    = [ -90,  90]
DEFAULT_RANGE_OFFSET = 0.00000000000000000000001
DEFAULT_AXIS         = [-180, 180, -90, 90]

def plot_mapped(data, baseMapInstance, title,
                vmin=None, vmax=None,
                boundingAxis=DEFAULT_AXIS,
                fillValue=DEFAULT_FILL_VALUE,
                colorMap=cm.jet) :
    
    # turn data into a masked array that excludes fill values and collapse 3D arrays
    temp_mask = None
    if not data is None :
        # TODO, there's probably a better way to show 3D data
        if len(data.shape) > 2 :
            # calculate a rough mean
            data = numpy.nansum(data, axis=0) / numpy.sum(numpy.isfinite(data), axis=0)
            #data = data[0]
        
        temp_mask = numpy.isnan(data) if fillValue is numpy.nan else data == fillValue
        data      = numpy.ma.masked_where(temp_mask, data)
    
    lat_data, lon_data = make_lat_lon_grids((data.shape[1], data.shape[0]))
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # figure the range for the color bar
    levelsToUse = None
    if not (data is None) :
        
        minVal = vmin if vmin is not None else data.min() - DEFAULT_RANGE_OFFSET
        maxVal = vmax if vmax is not None else data.max() + DEFAULT_RANGE_OFFSET
        
        levelsToUse = numpy.linspace(minVal, maxVal, DEFAULT_LEVELS_NUM)
    
    # draw our data placed on a map
    draw_basic_features(baseMapInstance, boundingAxis)
    x, y = baseMapInstance(lon_data, lat_data)
    p    = baseMapInstance.contourf(x, y, data, levelsToUse, cmap=colorMap)
    
    # set the title
    axes.set_title(title)
    
    # show a generic color bar
    if data is not None :
        cbar = plt.colorbar(format='%.3g')

def plot_binary(data, title,
                fill_value=DEFAULT_FILL_VALUE,
                vmin=None, vmax=None):
    
    # TODO, there's probably a better way to show 3D data
    raw_data   = numpy.nansum(data, axis=0) if len(data.shape) > 2 else data
    
    # mask the data based on the fill value
    temp_mask   = numpy.isnan(raw_data) if fill_value is numpy.nan else raw_data == fill_value
    masked_data = numpy.ma.masked_where(temp_mask, raw_data)
    
    # plot the figure
    figure = plt.figure()
    axes   = figure.add_subplot(111)
    plt.imshow(masked_data, vmin=vmin, vmax=vmax)
    axes.set_title(title)
    plt.bone()
    plt.colorbar()

def create_basemap (axis=DEFAULT_AXIS, projection='cyl', resolution='i') :
    """
    Create an instance of basemap using either the specified axis info or the
    specified lon and lat info to pick the viewing area.
    the format of the axis is [lon min, lon max, lat min, lat max]
    where the min and max are for the entire area that you wish to show.
    
    Note: There are known viewing area problems with conic projections that
    may cause "rectangular" data to be clipped.
    """
    
    # pull out the longitude/latitude info
    lon_left   = axis[0] 
    lat_bottom = axis[2] 
    lon_right  = axis[1] 
    lat_top    = axis[3] 
    lon_mid    = (lon_left + lon_right ) / 2.
    lat_mid    = (lat_top  + lat_bottom) / 2.
    
    # make our basemap
    m = None
    if projection is 'ortho' :
        # orthographic projections require this call
        m = Basemap(resolution=resolution, area_thresh=10000., projection=projection,
                    lat_0=lat_mid, lon_0=lon_mid)
    else :
        # most of the other projections use this call
        m = Basemap(llcrnrlon=lon_left,llcrnrlat=lat_bottom,urcrnrlon=lon_right,urcrnrlat=lat_top,
                    resolution=resolution, area_thresh=10000., projection=projection,
                    lat_1=lat_mid,lon_0=lon_mid)
    
    return m, axis

def draw_basic_features(baseMapInstance, axis) :
    """
    Draw the basic outlines of the earth's features.
    """
    # draw the basic physical and geopolitical features
    baseMapInstance.drawcoastlines()
    baseMapInstance.drawcountries()
    baseMapInstance.drawstates()
    baseMapInstance.drawmapboundary()
    
    # pull out the longitude/latitude info
    lon_left   = axis[0] 
    lat_bottom = axis[2] 
    lon_right  = axis[1] 
    lat_top    = axis[3] 
    
    # draw the parallels and meridians
    parallels = numpy.arange(-80.,90.,abs(lat_top - lat_bottom) / 4.0)
    baseMapInstance.drawparallels(parallels,labels=[1,0,0,1])
    meridians = numpy.arange(0., 360.,abs(lon_left - lon_right) / 4.0)
    baseMapInstance.drawmeridians(meridians,labels=[1,0,0,1])    

def load_fbf (file_name, in_dir='.') :
    """
    load raw data from a flat binary file
    """
    
    # get the data from the file
    var_workspace = Workspace.Workspace(dir=in_dir)
    fbf_attr_name = file_name.split(".")[0]
    raw_data      = var_workspace[fbf_attr_name][:]
    
    return raw_data, fbf_attr_name

def save_last_plot (fbf_attr_name, dpi_to_use=DEFAULT_DPI) :
    """
    save the last figure plotted with plt to an appropriately named file
    """
    
    plt.savefig(PLOT_SUFFIX + "." + fbf_attr_name + ".png", dpi=dpi_to_use)
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
    
    parser.add_argument('-m', '--mapped', dest="do_plot_mapped", default=False, action="store_true",
                        help="Specify that the data should be plotted on a lon/lat grid with background maps")
    
    args = parser.parse_args()
    
    workspace = '.'
    binary_files = args.binary_files
    if not args.binary_files and not args.pattern:
        args.pattern = DEFAULT_FILE_PATTERN
    if args.pattern:
        workspace = os.path.split(args.pattern)[0]
        binary_files = [ os.path.split(x)[1] for x in glob(args.pattern) ]
    
    basemap_temp = None
    if args.do_plot_mapped :
        basemap_temp, _ = create_basemap()
    
    for bf in binary_files:
        print "Plotting '%s'" % (bf,)
        try:
            raw_data, var_name = load_fbf (bf)
            
            if args.do_plot_mapped :
                plot_mapped(raw_data, basemap_temp, var_name + " data",
                            vmin=args.vmin, vmax=args.vmax)
            else :
                plot_binary(raw_data, var_name + " data",
                            fill_value=args.fill_value,
                            vmin=args.vmin, vmax=args.vmax)
            
            save_last_plot (var_name, dpi_to_use=args.dpi)
            
        except StandardError as e:
            print "Could not plot '%s'" % (bf,)
            if hasattr(e, "msg"): print e,e.msg
            else: print e

if __name__ == "__main__":
    import sys

    sys.exit(main())

