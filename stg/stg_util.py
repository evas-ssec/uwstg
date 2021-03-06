#!/usr/bin/env python
# encoding: utf-8
"""
Handle useful support functions for stg.

:author:       Eva Schiffer (evas)
:contact:      evas@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
:date:         Feb 2014
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

from constants import *

import numpy
import os, logging

LOG = logging.getLogger(__name__)

def make_lat_lon_grids ((lat_size, lon_size), lat_min=-90, lat_max=90, lon_min=-180, lon_max=180) :
    
    # build our lon/lat from the data shape
    lon_row  = numpy.linspace(lon_min, lon_max, lon_size)
    lon_data = numpy.tile(lon_row, (lat_size, 1))
    #lon_data = numpy.transpose(lon_data)
    
    lat_row  = numpy.linspace(lat_min, lat_max, lat_size)
    lat_data = numpy.tile(lat_row, (lon_size, 1))
    lat_data = numpy.transpose(lat_data)
    
    return lat_data, lon_data

def make_index_grid ((lat_size, lon_size)) :
    
    temp_lat, temp_lon = make_lat_lon_grids((lat_size, lon_size), lat_min=0, lat_max=lat_size-1, lon_min=0, lon_max=lon_size-1)
    temp_lat = temp_lat.astype(numpy.int)
    temp_lon = temp_lon.astype(numpy.int)
    
    return temp_lat, temp_lon

def clean_path(string_path) :
    """
    Return a clean form of the path without any '.', '..', or '~'
    """
    clean_path = None
    if string_path is not None :
        clean_path = os.path.abspath(os.path.expanduser(string_path))

    return clean_path

def setup_dir_if_needed(dir_path, description_name) :
    """
    create the directory if that is needed, if not don't
    """
    if not (os.path.isdir(dir_path)) :
        LOG.info("Specified " + description_name + " directory (" + dir_path + ") does not exist.")
        LOG.info("Creating " + description_name + " directory.")
        os.makedirs(dir_path)
