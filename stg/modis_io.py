#!/usr/bin/env python
# encoding: utf-8
"""
Handle parsing input from Modis files.

:author:       Eva Schiffer (evas)
:contact:      evas@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
:date:         Jan 2014
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

from constants import *

import sys
import logging

import numpy
from pyhdf.SD import SD,SDC, SDS, HDF4Error

from stg.constants import *
import stg.modis_guidebook as modis_guidebook

LOG = logging.getLogger(__name__)

# the line between day and night for our day/night masks (in solar zenith angle degrees)
DAY_NIGHT_LINE_DEGREES = 84.0

def open_file (file_path) :
    """
    given a file path that is a modis file, open it
    """
    
    file_object = SD(file_path, SDC.READ)
    
    return file_object

def close_file (file_object) :
    """
    given a file object, close it
    """
    
    file_object.end()

def load_aux_data (file_path, minimum_scan_angle, file_object=None) :
    """
    load the auxillary data and process the appropriate masks from it
    """
    
    # make our return structure
    aux_data_sets = { }
    
    # load the longitude and latitude
    file_object, aux_data_sets[LON_KEY] = load_variable_from_file (modis_guidebook.LONGITUDE_NAME,
                                                                   file_path=file_path, file_object=file_object)
    file_object, aux_data_sets[LAT_KEY] = load_variable_from_file (modis_guidebook.LATITUDE_NAME,
                                                                   file_path=file_path, file_object=file_object)
    
    # load the angles to make masks
    file_object, solar_zenith_data_temp = load_variable_from_file (modis_guidebook.SOLAR_ZENITH_NAME,
                                                                   file_path=file_path, file_object=file_object)
    file_object, sat_zenith_data_temp   = load_variable_from_file (modis_guidebook.SENSOR_ZENITH_NAME,
                                                                   file_path=file_path, file_object=file_object)
    
    # transform the satellite zenith to scan angle
    scan_angle_data_temp = satellite_zenith_angle_to_scan_angle(sat_zenith_data_temp)
    
    # build the day and night masks
    ok_scan_angle                 =  scan_angle_data_temp   <= minimum_scan_angle
    aux_data_sets[DAY_MASK_KEY]   = (solar_zenith_data_temp <  DAY_NIGHT_LINE_DEGREES) & ok_scan_angle
    aux_data_sets[NIGHT_MASK_KEY] = (solar_zenith_data_temp >= DAY_NIGHT_LINE_DEGREES) & ok_scan_angle
    
    return file_object, aux_data_sets

# FUTURE, the data type needs to be handled differently
def load_variable_from_file (variable_name, file_path=None, file_object=None,
                             fill_value_name=modis_guidebook.FILL_VALUE_ATTR_NAME,
                             scale_name=modis_guidebook.SCALE_ATTR_NAME,
                             offset_name=modis_guidebook.ADD_OFFSET_ATTR_NAME,
                             data_type_for_output=numpy.float32) :
    """
    load a given variable from a file path or file object
    """
    
    if file_path is None and file_object is None :
        raise ValueError("File path or file object must be given to load file.")
    if file_object is None :
        file_object = open_file(file_path)
    
    variable_names = file_object.datasets().keys()
    if variable_name not in variable_names :
        raise ValueError("Variable " + str(variable_name) +
                         " is not present in file " + str(file_path) + " .")
    
    # defaults
    scale_factor = 1.0
    add_offset = 0.0
    data_type = None
    
    # get the variable object and use it to
    # get our raw data and scaling info
    variable_object = file_object.select(variable_name)
    raw_data_copy   = variable_object[:]
    raw_data_copy   = raw_data_copy.astype(data_type_for_output) if data_type_for_output is not None else raw_data_copy
    temp_attrs      = variable_object.attributes()
    try :
        scale_factor, scale_factor_error, add_offset, add_offset_error, data_type = SDS.getcal(variable_object)
    except HDF4Error:
        # load just the scale factor and add offset information by hand
        if offset_name in temp_attrs.keys() :
            add_offset = temp_attrs[offset_name]
            data_type = numpy.dtype(type(add_offset))
        if scale_name in temp_attrs.keys() :
            scale_factor = temp_attrs[scale_name]
            data_type = numpy.dtype(type(scale_factor))
    
    # get the fill value
    fill_value = temp_attrs[fill_value_name] if fill_value_name in temp_attrs.keys() else numpy.nan
    # change the fill value to numpy.nan
    fill_mask = raw_data_copy == fill_value
    if fill_value is not numpy.nan :
        raw_data_copy[fill_mask] = numpy.nan
        fill_value = numpy.nan
    
    # we got all the info we need about that file
    SDS.endaccess(variable_object)
    
    # don't do lots of work if we don't need to scale things
    if (scale_factor == 1.0) and (add_offset == 0.0) :
        return file_object, raw_data_copy
    
    # if we don't have a data type something strange has gone wrong
    assert(data_type is not None)
    
    # create the scaled version of the data
    scaled_data_copy                = raw_data_copy.copy()
    scaled_data_copy = unscale_data(scaled_data_copy, fill_mask=fill_mask,
                                    scale_factor=scale_factor, offset=add_offset)
    
    return file_object, scaled_data_copy 

def unscale_data (data, fill_mask=None, scale_factor=None, offset=None) :
    """unscale the given data
    
    data is modified in place and fill values will not be changed
    if a scale factor or offset is given as None (or not given) it will not be applied
    
    the general formula for unscaling data is:
    
            final_data = scale_factor * (input_data - offset)
    
    """
    
    to_return = data
    
    not_fill_mask = ~fill_mask
    
    # if we have an offset use it to offset the data
    if (offset is       not None) and (offset       is not 0.0) :
        to_return[not_fill_mask] -= offset
    
    # if we found a scale use it to scale the data
    if (scale_factor is not None) and (scale_factor is not 1.0) :
        to_return[not_fill_mask] *= scale_factor
    
    return to_return

# FUTURE, will this be used for other satellites? should it move up to the io_manager?
def satellite_zenith_angle_to_scan_angle (sat_zenith_data) :
    """
    given a set of satellite zenith angles, calculate the equivalent scan angles
    
    Note: This comes directly from Nadia's satz2scang function.
    """
    
    # some constants
    re = 6371.03;
    altkm = 825;
    fac = re / (re + altkm);
    dtr = 0.01745329;
    
    # do the angle calculations
    arg_data        = sat_zenith_data * dtr
    ang_data        = numpy.sin(numpy.sin(arg_data) * fac)
    scan_angle_data = ang_data / dtr
    
    return scan_angle_data

def main():
    import optparse
    from pprint import pprint
    usage = """
%prog [options] filename1.hdf

"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
            help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_option('-r', '--no-read', dest='read_hdf', action='store_false', default=True,
            help="don't read or look for the hdf file, only analyze the filename")
    (options, args) = parser.parse_args()
    
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[min(3, options.verbosity)])
    
    LOG.info("Currently no command line tests are set up for this module.")

if __name__ == '__main__':
    sys.exit(main())
