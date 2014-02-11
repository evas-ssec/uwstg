#!/usr/bin/env python
# encoding: utf-8
"""
Handle parsing input from Rich Frey's binary CTP files.

:author:       Nick Bearson (nickb)
:contact:      nickb@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
:date:         Feb 2014
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

from constants import *

import sys
import logging

import numpy

from stg.constants import *
import stg.ctp_guidebook as ctp_guidebook

LOG = logging.getLogger(__name__)

# the line between day and night for our day/night masks (in solar zenith angle degrees)
DAY_NIGHT_LINE_DEGREES = 84.0

def open_file (file_path) :
    """
    given a file path that is a modis file, open it
    """
    ctp_names = [
             ctp_guidebook.LATITUDE_NAME              ,
             ctp_guidebook.LONGITUDE_NAME             ,
             ctp_guidebook.CLOUD_TOP_PRESS_NAME       ,
             ctp_guidebook.CLOUD_TOP_HEIGHT_NAME      ,
             ctp_guidebook.CLOUD_TOP_TEMP_NAME        ,
             ctp_guidebook.EFFECTIVE_CLOUD_AMOUNT_NAME,
             ctp_guidebook.METHOD_FLAG_NAME           ,
             ctp_guidebook.DAY_NIGHT_FLAG_NAME        ,
             ctp_guidebook.DIRECTION_FLAG_NAME        ,
             ctp_guidebook.VIEWING_ZENITH_NAME        ,
             ctp_guidebook.SCAN_LINE_TIME_NAME        ,
             ctp_guidebook.CLOUD_FRACTION_NAME        ,
             ctp_guidebook.LAND_FRACTION_NAME         ,
             ctp_guidebook.RESULTS_FLAG_NAME          ,
             ctp_guidebook.UTLS_FLAG_NAME             ,
    ]

    ctp_formats = ['(1100,56)f4'] * 15   # they're all the same

    ctp_type = numpy.dtype({'names' : ctp_names,
                            'formats' : ctp_formats})

    file_object = numpy.fromfile(file_path, dtype=ctp_type) # FIXME: will we have a problem from treating a numpy array as a file object? passing it around?

    return file_object

def close_file (file_object) :
    """
    given a file object, close it
    """
    del file_object


def load_aux_data (file_path, minimum_scan_angle, file_object=None) :
    """
    load the auxillary data and process the appropriate masks from it
    """
    
    # make our return structure
    aux_data_sets = { }
    
    # load the longitude and latitude
    file_object, aux_data_sets[LON_KEY] = load_variable_from_file (ctp_guidebook.LONGITUDE_NAME,
                                                                   file_path=file_path, file_object=file_object)
    file_object, aux_data_sets[LAT_KEY] = load_variable_from_file (ctp_guidebook.LATITUDE_NAME,
                                                                   file_path=file_path, file_object=file_object)
    
    # load the day/night flag to make day/night mask
    file_object, day_night_flag         = load_variable_from_file (ctp_guidebook.DAY_NIGHT_FLAG_NAME,
                                                                   file_path=file_path, file_object=file_object)

    # build the day and night masks
    aux_data_sets[DAY_MASK_KEY]   = (day_night_flag == 1)
    aux_data_sets[NIGHT_MASK_KEY] = (day_night_flag == 2)
    
    return file_object, aux_data_sets

# FUTURE, the data type needs to be handled differently
def load_variable_from_file (variable_name, file_path=None, file_object=None,
                             fill_value_name=None,
                             scale_name=None,
                             offset_name=None,
                             data_type_for_output=numpy.float32) :
    """
    load a given variable from a file path or file object
    """
    if file_path is None and file_object is None :
        raise ValueError("File path or file object must be given to load file.")
    if file_object is None :
        file_object = open_file(file_path)

    data = file_object[variable_name]

# FIXME: everything appears to have negative fill values (aside from lat/lon), verify with RichF
    if variable_name not in [ctp_guidebook.LONGITUDE_NAME, ctp_guidebook.LATITUDE_NAME]:
      fill_mask = (data < 0)
      data[fill_mask] = numpy.nan

    data_to_return = data.astype(data_type_for_output) if data_type_for_output is not None else data
    return file_object, data


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
