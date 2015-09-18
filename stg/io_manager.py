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

import sys, logging, os, re

import numpy

import stg.modis_guidebook   as modis_guidebook
import stg.modis_io          as modis_io

import stg.ctp_guidebook as ctp_guidebook
import stg.ctp_io        as ctp_io

import keoni.fbf as fbf

LOG = logging.getLogger(__name__)

# these are suffixes used for temporary files while space gridding
EXPECTED_TEMP_SUFFIXES    = [
                             "%s"   + "-" + TEMP_SUFFIX_KEY,
                             "%s"   + "-" + DENSITY_SUFFIX + "-" + TEMP_SUFFIX_KEY,
                             "%s"   + "-" + NOBS_SUFFIX + "-" + TEMP_SUFFIX_KEY,
                            ]
# these are suffixes used for the final, packed files from space gridding
EXPECTED_SPACE_OUT_SUFFIXES   = [
                                 "%s"   + "-" + DAILY_SPACE_SUFFIX_KEY,
                                 "%s"   + "-" + NOBS_SUFFIX + "-" + DAILY_SPACE_SUFFIX_KEY,
                                ]

# the strftime format for date stamping our files
DATE_STAMP_FORMAT         = "%Y%m%d"

def get_list_of_suffixes (gridding_type, suffix_type) :
    """make a list of suffixes based on the possible times and the type of files these are

    gridding_type is expected to be: DAILY_SPACE_TYPE, DAILY_TIME_TYPE, MULTIDAY_TIME_TYPE, or NOBS_LUT_TYPE

    suffix_type is expected to be: TEMP_FILE_TYPE, FINAL_OUT_FILE_TYPE, or ALL_FILES_TYPE
    """

    to_return = [ ]

    if gridding_type == DAILY_SPACE_TYPE :

        # add the temp suffixes to the list if desired
        if suffix_type == TEMP_FILE_TYPE or suffix_type == ALL_FILES_TYPE :

            for time_key in TIME_SETS :
                for suffix in EXPECTED_TEMP_SUFFIXES :
                    to_return.append(suffix % time_key)

        if suffix_type == FINAL_OUT_FILE_TYPE or suffix_type == ALL_FILES_TYPE :

            for time_key in TIME_SETS :
                for suffix in EXPECTED_SPACE_OUT_SUFFIXES :
                    to_return.append(suffix % time_key)

    elif gridding_type == DAILY_TIME_TYPE :
        LOG.debug("Daily time gridded file suffixes are not yet provided by this method.") # TODO
    elif gridding_type == MULTIDAY_TIME_TYPE :
        LOG.debug("Multi-day time gridded file suffixes are not yet provided by this method.") # TODO
    elif gridding_type == NOBS_LUT_TYPE :
        LOG.debug("Num-obs LUT file suffixes are not yet provided by this method.") # TODO
    else :
        LOG.debug("Unexpected file type (" + gridding_type + ") requested in get_list_of_suffixes.")

    return to_return

def open_file (file_path) :
    """
    given a file path, open it
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

def get_expected_abstract_sets (satellite_constant, separate_day_night=True) :
    
    expected_data_sets = { }
    
    if satellite_constant == INST_MODIS :
        expected_data_sets = modis_io.get_abstract_data_sets (do_separate_day_night=separate_day_night)
    # FUTURE, needs a statement for ctp
    
    return expected_data_sets

def get_expected_data_sets_from_aux_data (satellite_constant, aux_data, do_separate_day_night=True) :
    """given aux data in the form returned by load_aux_data and the grid degrees constant, return the data sets to be processed
    
    Each data set is defined by a constant name, a mask to select that set, it's expected suffixes for temporary density/nobs/data
    and it's expected suffix for the final output data/nobs
    """
    
    expected_data_sets = { }
    
    if satellite_constant == INST_MODIS :
        expected_data_sets = modis_io.determine_data_sets(aux_data, do_separate_day_night=do_separate_day_night)
    # FUTURE, needs a statement for ctp
    
    return expected_data_sets

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
    temp_path = os.path.join(os.path.abspath(os.path.expanduser(output_path)), temp_file)
    temp_file_obj = open(temp_path, file_permissions)
    data_array.astype(data_type).tofile(temp_file_obj)
    temp_file_obj.close()

def build_name_stem (variable_name, date_time=None, satellite=None, suffix=None) :
    """given information on what's in the file, build a file stem
    if there's extra info like the date time, satellite, algorithm name, or a suffix
    include that in the file stem as well
    
    the name format is:
            datestamp_satellite_variablename_suffix
    """
    
    # the basic stem name is just the variable
    stem_name = variable_name
    
    # if we have a satellite, add that to the beginning
    stem_name = satellite + "_" + stem_name if satellite is not None else stem_name

    # if we have a date time, add a time stamp at the beginning
    stem_name = date_time.strftime(DATE_STAMP_FORMAT) + "_" + stem_name if date_time is not None else stem_name

    # if we have a suffix, add that too
    stem_name = stem_name + "_" + suffix if suffix is not None else stem_name
    
    return stem_name
    
    # date_stamp + "_" + var_name + suffix

def parse_flatfile_name (flatfile_name) :
    """given a flat file name, parse out the things we expect from the stem

    The shape of the variable will also be returned as a convenience.

    Note: there is some ambiguity in the flat file names because some variables use "_", so
    this method is only guaranteed to work correctly when the date_time, satellite, and suffix were
    included in the original flat file name stem
    """

    LOG.debug("Parsing flat file name: " + flatfile_name)

    var_name    = None
    datetimestr = None
    satellite   = None
    suffix      = None
    var_shape   = None
    file_type   = None

    # strip off the shape and data type from the end
    temp_name = flatfile_name.split(".")
    if len(temp_name) >= 4 :
        #data_type = temp_name[1] # TODO, would it improve the algorithm to use this?
        var_shape = tuple(reversed(temp_name[2:]))
    temp_name = temp_name[0] # extract the name stem

    # parse the part of the name not related to the flat file structure
    _split = temp_name.split("_")
    if len(_split) >= 4 :

        # detect if there is a date on the front, if so, strip that
        try :
            int(_split[0])
            datetimestr = _split[0]
            _split = _split[1:]
        except ValueError :
            LOG.debug("Unable to strip date from flat file name: " + flatfile_name)

        # detect if there is a satellite on the front, if so, strip that
        temp_sat = _split[0]
        if temp_sat in ALL_SATS :
            satellite = temp_sat
            _split = _split[1:]
        else :
            LOG.debug("Unable to identify satellite name from flat file name: " + flatfile_name)

        # detect if this is a space gridded file
        if _split[-1].find(DAILY_SPACE_SUFFIX_KEY) >= 0 :
            file_type = DAILY_SPACE_TYPE
            suffix = _split[-1]
            _split = _split[:-1]
        # detect if this is a daily time gridded file
        elif _split[-1].find(DAILY_TIME_SUFFIX_KEY) >= 0 :
            file_type = DAILY_TIME_TYPE
            suffix = _split[-1]
            _split = _split[:-1]
        # detect if this is a multi-day time gridded file
        elif _split[-1].find(MULTI_TIME_SUFFIX_KEY) >= 0 :
            file_type = MULTIDAY_TIME_TYPE
            suffix = _split[-1]
            _split = _split[:-1]
        # detect if this is a nobs look up table
        elif _split[-1].find(NOBS_LUT_SUFFIX) >= 0 :
            file_type = NOBS_LUT_TYPE
            suffix = _split[-1]
            _split = _split[:-1]
        else :
            LOG.debug("Unable to determine file type from suffix in flat file name: " + flatfile_name)

        # the remaining part of the name should be the variable name
        if len(_split) > 0 :
            var_name = ""
            for piece in _split :
                var_name = var_name + "_" + piece
            var_name = var_name[1:] # take off the first, unneeded underscore

    return file_type, var_name, datetimestr, satellite, suffix, var_shape

def get_date_stamp_from_file_name (file_name) :
    """given a file name starting with a name stem created by build_name_stem, determine the date stamp for the file
    """
    
    date_stamp_to_return = None
    
    # make sure we only have the stem
    stem = file_name.split('.')[0]
    
    # break the stem up by underscores
    split_stem = stem.split('_')
    for section in split_stem :
        
        if re.match(r'\d\d\d\d\d\d\d\d', section) :
            
            date_stamp_to_return = section
    
    return date_stamp_to_return

def organize_space_gridded_files (file_name_list) :
    """organize the files into sets based on their names
    
    Note: the type of the files is determined by the first file, so the
    function will not handle mixed sets of files properly
    """
    
    to_return = { }
    
    first_file = file_name_list[0]
    if modis_guidebook.is_MODIS_flat_file(first_file) :
        to_return = modis_io.organize_space_gridded_files(file_name_list)
    # FUTURE, needs a statement for ctp
    
    return to_return

def main():
    import optparse
    #from pprint import pprint
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
