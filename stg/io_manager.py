#!/usr/bin/env python
# encoding: utf-8
"""
Handle parsing input from files.

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
import os

import numpy

import stg.general_guidebook as general_guidebook

import stg.modis_guidebook   as modis_guidebook
import stg.modis_io          as modis_io

import stg.ctp_guidebook as ctp_guidebook
import stg.ctp_io        as ctp_io

import keoni.fbf as fbf

LOG = logging.getLogger(__name__)

# these are suffixes used for temporary files
DAY_TEMP_SUFFIX           = "_daytemp"
NIGHT_TEMP_SUFFIX         = "_nighttemp"
DAY_DENSITY_TEMP_SUFFIX   = "_daydensitytemp"
NIGHT_DENSITY_TEMP_SUFFIX = "_nightdensitytemp"
DAY_NOBS_TEMP_SUFFIX      = "_daynobstemp"
NIGHT_NOBS_TEMP_SUFFIX    = "_nightnobstemp"
EXPECTED_TEMP_SUFFIXES    = [DAY_TEMP_SUFFIX,         NIGHT_TEMP_SUFFIX,
                             DAY_DENSITY_TEMP_SUFFIX, NIGHT_DENSITY_TEMP_SUFFIX,
                             DAY_NOBS_TEMP_SUFFIX,    NIGHT_NOBS_TEMP_SUFFIX]

# these are suffixes used for the final, packed files
DAY_SUFFIX                = "_dayfinal"
NIGHT_SUFFIX              = "_nightfinal"
DAY_NOBS_SUFFIX           = "_daynobsfinal"
NIGHT_NOBS_SUFFIX         = "_nightnobsfinal"
EXPECTED_FINAL_SUFFIXES   = [DAY_SUFFIX,      NIGHT_SUFFIX,
                             DAY_NOBS_SUFFIX, NIGHT_NOBS_SUFFIX]

# all the suffixes we can produce
ALL_EXPECTED_SUFFIXES     = [DAY_TEMP_SUFFIX,         NIGHT_TEMP_SUFFIX,
                             DAY_DENSITY_TEMP_SUFFIX, NIGHT_DENSITY_TEMP_SUFFIX,
                             DAY_NOBS_TEMP_SUFFIX,    NIGHT_NOBS_TEMP_SUFFIX,
                             DAY_SUFFIX,              NIGHT_SUFFIX,
                             DAY_NOBS_SUFFIX,         NIGHT_NOBS_SUFFIX]

# the strftime format for date stamping our files
DATE_STAMP_FORMAT         = "%Y%m%d"

def open_file (file_path) :
    """
    given a file path that is a modis file, open it
    """
    
    file_object = None
    
    if modis_guidebook.is_MODIS_file(file_path) :
        file_object = modis_io.open_file(file_path)

    if ctp_guidebook.is_CTP_file(file_path):
        file_object = ctp_io.open_file(file_path)

    return file_object

def close_file (file_path, file_object) :
    """
    given a file object, close it
    """
    
    if modis_guidebook.is_MODIS_file(file_path) :
        modis_io.close_file(file_object)

    if ctp_guidebook.is_CTP_file(file_path):
        ctp_io.close_file(file_object)

def load_aux_data (file_path, minimum_scan_angle, file_object=None) :
    """
    load the auxillary data for a given file
    """
    
    temp_aux_data = None
    if modis_guidebook.is_MODIS_file(file_path) :
        file_object, temp_aux_data = modis_io.load_aux_data(file_path,
                                                            minimum_scan_angle,
                                                            file_object=file_object)

    if ctp_guidebook.is_CTP_file(file_path):
        file_object, temp_aux_data = ctp_io.load_aux_data(file_path,
                                                            minimum_scan_angle,
                                                            file_object=file_object)

    return file_object, temp_aux_data

def load_variable_from_file (variable_name, file_path=None, file_object=None,
                             data_type_for_output=numpy.float32) :
    """
    load a given variable from a file path or file object
    """
    
    temp_data = None
    
    if modis_guidebook.is_MODIS_file(file_path) :
        file_object, temp_data = modis_io.load_variable_from_file (variable_name,
                                                                   file_path=file_path,
                                                                   file_object=file_object,
                                                                   data_type_for_output=data_type_for_output)
    
    if ctp_guidebook.is_CTP_file(file_path) :
        file_object, temp_data = ctp_io.load_variable_from_file (variable_name,
                                                                   file_path=file_path,
                                                                   file_object=file_object,
                                                                   data_type_for_output=data_type_for_output)

    return file_object, temp_data

def save_data_to_file (stem_name, grid_shape, output_path, data_array, data_type, file_permissions="a") :
    """
    save numpy data to the appropriate file name given the information about the stem and path
    
    by default this will append to the end of a file, you can instead pass in a different set of
    file_permissions, using the standard permissions from python's open function
    """
    
    temp_file = fbf.filename(stem_name, data_type, shape=grid_shape)
    temp_path = os.path.join(output_path, temp_file)
    temp_file_obj = open(temp_path, 'a')
    data_array.astype(data_type).tofile(temp_file_obj)
    temp_file_obj.close()

def build_name_stem (variable_name, date_time=None, satellite=None, algorithm=None, suffix=None) :
    """given information on what's in the file, build a file stem
    if there's extra info like the date time, satellite, algorithm name, or a suffix
    include that in the file stem as well
    
    the name format is:
            satellite_algorithm_datestamp_variablename_suffix
    """
    
    # the basic stem name is just the variable
    stem_name = variable_name
    
    # if we have a date time, add a time stamp at the beginning
    stem_name = date_time.strftime(DATE_STAMP_FORMAT) + "_" + stem_name if date_time is not None else stem_name
    
    # if we have an algorithm prefix add that
    stem_name = algorithm + "_" + stem_name if algorithm is not None else stem_name
    
    # if we have a satellite, add that to the beginning
    stem_name = satellite + "_" + stem_name if satellite is not None else stem_name
    
    # if we have a suffix, add that too
    stem_name = stem_name + suffix if suffix is not None else stem_name
    
    return stem_name
    
    # date_stamp + "_" + var_name + suffix

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
