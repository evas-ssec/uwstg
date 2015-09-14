#!/usr/bin/env python
# encoding: utf-8
"""
This module handles the basic space gridding algorithm.

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

import numpy
import logging

from stg.constants import *
from stg.stg_util  import make_index_grid

LOG = logging.getLogger(__name__)

# this represents the time range of data that we will consider the "same overpass"
# for the purposes of putting it in the same grid cell when space gridding data
SAME_TIME_RANGE_SECONDS = 1800.0

def calculate_index_from_nav_data (lat_data, lon_data, grid_degrees) :
    """given the lon/lat data, calculate where the elements will be space gridded

    note, grid_degrees defines how large a grid cell will be, ie. the coarseness of the grid
    """
    
    # calculate the indexes for lon and lat
    lon_index = numpy.round((lon_data + 180.0) / grid_degrees) % (360.0 / grid_degrees)
    lat_index = numpy.round((lat_data +  90.0) / grid_degrees) % (180.0 / grid_degrees)
    
    return lat_index, lon_index

def space_grid_data (grid_lat_size, grid_lon_size, data, lat_indexes, lon_indexes ) :
    """given lon/lat indexes, data, and the grid size, sort the data into a space grid
    
    returns the filled space grid (empty space is NaN values), a density map of where the data is, and the size of the deepest bucket
    """
    
    #if data.size > 0 :
    #    print ("data range: " + str(data.min()) + " " + str(data.max()))
    
    space_grid_shape = (grid_lat_size, grid_lon_size) # I've confirmed with Nadia that this is the correct order
    
    # create the density map and figure out how dense the data will be
    # FUTURE, figure out how to do this in native numpy ops instead of loops
    density_map = numpy.zeros(space_grid_shape)
    nobs_map    = numpy.zeros(space_grid_shape)
    for index in range(data.size) :
        nobs_map[lat_indexes[index], lon_indexes[index]] += 1
        if numpy.isfinite(data[index]) :
            density_map[lat_indexes[index], lon_indexes[index]] += 1
    max_depth = numpy.max(density_map)
    
    #print ("max depth: " + str(max_depth))
    
    # create the space grids for this variable
    space_grid = numpy.ones((max_depth, grid_lat_size, grid_lon_size), dtype=numpy.float32) * numpy.nan #TODO, dtype should be set dynamically
    temp_depth = numpy.zeros(space_grid_shape)
    
    # put the variable data into the space grid
    # FUTURE, figure out how to do this in native numpy ops instead of loops
    for index in range(data.size) :
        if numpy.isfinite(data[index]) :
            depth = temp_depth[lat_indexes[index], lon_indexes[index]]
            space_grid[depth,  lat_indexes[index], lon_indexes[index]] = data[index]
            temp_depth[        lat_indexes[index], lon_indexes[index]] += 1

    return space_grid, density_map, nobs_map, max_depth

def pack_space_grid (data_array, density_array) :
    """
    given a sparse array of space gridded data and a density array
    where every slice represents the density of one "layer" in the
    sparse array, create a well packed version of the sparse array
    and return it
    """
    
    # figure out the maximum depth of the well packed data based on the density map
    max_depth = numpy.max(numpy.sum(density_array, axis=0))
    
    # create the final data array at the right depth
    final_data = numpy.ones((max_depth + 1, data_array.shape[1], data_array.shape[2]), dtype=data_array.dtype) * numpy.nan
    
    LOG.debug("  original data shape: " + str(data_array.shape))
    LOG.debug("  final data shape:    " + str(final_data.shape))
    
    # lat_data, lon_data = make_index_grid(data_array.shape[1:3])
    
    temp_index = numpy.zeros((data_array.shape[1], data_array.shape[2]), dtype=numpy.int)
    for depth in range(data_array.shape[0]) :
        
        this_slice = data_array[depth]
        valid_mask = numpy.isfinite(this_slice)
        
        # final_data[temp_index, lat_data, lon_data] = this_slice

        # FUTURE, can we do this directly with numpy instead of looping?
        for row in range(this_slice.shape[0]) :
            for col in range(this_slice.shape[1]) :
                if valid_mask[row, col] :
                    final_data[temp_index[row, col], row, col] = this_slice[row, col]
        
        
        temp_index += valid_mask
    
    return final_data[0:-1]


