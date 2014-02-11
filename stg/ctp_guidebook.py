#!/usr/bin/env python
# encoding: utf-8
"""
Provide information about products the system expects in Rich Frey's per-orbit binary CTP files, how
they're stored in the files, and how to manage them during the process of regridding.

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
import os
import logging
import numpy
from datetime import datetime

LOG = logging.getLogger(__name__)

# variable names expected in the files

LATITUDE_NAME               = 'Latitude'
LONGITUDE_NAME              = 'Longitude'
CLOUD_TOP_PRESS_NAME        = 'Cloud_Top_Pressure'
CLOUD_TOP_HEIGHT_NAME       = 'Cloud_Top_Height'
CLOUD_TOP_TEMP_NAME         = 'Cloud_Top_Temperature'
EFFECTIVE_CLOUD_AMOUNT_NAME = 'Effective_Cloud_Amount'
METHOD_FLAG_NAME            = 'Retrieval_Method_Flag'
DAY_NIGHT_FLAG_NAME         = 'Day_Night_Flag'
DIRECTION_FLAG_NAME         = 'Direction_Flag'
VIEWING_ZENITH_NAME         = 'Viewing_Zenith'
SCAN_LINE_TIME_NAME         = 'Scan_Line_Time'
CLOUD_FRACTION_NAME         = 'Cloud_Fraction'
LAND_FRACTION_NAME          = 'Land_Fraction'
RESULTS_FLAG_NAME           = 'Results_Flag'
UTLS_FLAG_NAME              = 'UTLS_Flag'

# a map of what the caller calls variables to variable names in files
CALLER_VARIABLE_MAP     = {
                            'latitude':                     LATITUDE_NAME,
                            'lat':                          LATITUDE_NAME,

                            'longitude':                    LONGITUDE_NAME,
                            'lon':                          LONGITUDE_NAME,

                            'cloud top pressure':           CLOUD_TOP_PRESS_NAME,
                            'pressure':                     CLOUD_TOP_PRESS_NAME,

                            'cloud top height':             CLOUD_TOP_HEIGHT_NAME,
                            'height':                       CLOUD_TOP_HEIGHT_NAME,


                            'cloud top temperature':        CLOUD_TOP_TEMP_NAME,
                            'ctt5km':                       CLOUD_TOP_TEMP_NAME,
                            'cloud_top_temperature':        CLOUD_TOP_TEMP_NAME,

                            'amount':                       EFFECTIVE_CLOUD_AMOUNT_NAME,
                            'cloud amount':                 EFFECTIVE_CLOUD_AMOUNT_NAME,
                            'effective cloud amount':       EFFECTIVE_CLOUD_AMOUNT_NAME,

                          }


# just do everything by default 
EXPECTED_VARIABLES_IN_FILE = set([LATITUDE_NAME              ,
                                  LONGITUDE_NAME             ,
                                  CLOUD_TOP_PRESS_NAME       ,
                                  CLOUD_TOP_HEIGHT_NAME      ,
                                  CLOUD_TOP_TEMP_NAME        , 
                                  EFFECTIVE_CLOUD_AMOUNT_NAME, 
                                  METHOD_FLAG_NAME           , 
                                  DAY_NIGHT_FLAG_NAME        , 
                                  DIRECTION_FLAG_NAME        , 
                                  VIEWING_ZENITH_NAME        , 
                                  SCAN_LINE_TIME_NAME        , 
                                  CLOUD_FRACTION_NAME        , 
                                  LAND_FRACTION_NAME         , 
                                  RESULTS_FLAG_NAME          , 
                                  UTLS_FLAG_NAME             ,]) 
 

# TODO, move this up to the general_guidebook
def _clean_off_path_if_needed(file_name_string) :
    """
    remove the path from the file if necessary
    """
    
    return os.path.basename(file_name_string)

def is_CTP_file (file_name_string) :
    """determine if a file name is the right pattern to represent a CTP file
    if the file_name_string matches how we expect CTP files to look return
    TRUE else will return FALSE
    """
    
    temp_name_string = _clean_off_path_if_needed(file_name_string)
    
    return (temp_name_string.endswith('ctp.bin'))

def parse_datetime_from_filename (file_name_string) :
    """parse the given file_name_string and create an appropriate datetime object
    that represents the datetime indicated by the file name; if the file name does
    not represent a pattern that is understood, None will be returned
    """
    
    temp_name_string = _clean_off_path_if_needed(file_name_string)

    temp = temp_name_string.split('.')
    datetime_to_return = datetime.strptime(temp[3] + temp[4], 'D%y%jS%H%M')

    return datetime_to_return

def get_satellite_from_filename (data_file_name_string) :
    """given a file name, figure out which satellite it's from
    if the file does not represent a known satellite name
    configuration None will be returned
    """
    
    temp_name_string = _clean_off_path_if_needed(data_file_name_string)
    temp = temp_name_string.split('.')

    SAT_LOOKUP = {'NK' : SAT_NOAA_15,
                  'NL' : SAT_NOAA_16,
                  'NM' : SAT_NOAA_17,
                  'NN' : SAT_NOAA_18,
                  'NP' : SAT_NOAA_19,
                  'NK' : SAT_NOAA_15,
                  'M2' : SAT_METOP_A, # yes, M2 really does = metop-a
                  'M1' : SAT_METOP_B,
    }

    satellite_to_return = SAT_LOOKUP[temp[3]]
    instrument_to_return = INST_HIRS

    return satellite_to_return, instrument_to_return

def get_variable_names (user_requested_names) :
    """get a list of variable names we expect to process from the file
    """
    
    var_names = set( )
    
    if len(user_requested_names) <= 0 :
        var_names.update(EXPECTED_VARIABLES_IN_FILE)
    else :

        for user_name in user_requested_names :
            if user_name in CALLER_VARIABLE_MAP.keys() :
                var_names.update(set([CALLER_VARIABLE_MAP[user_name]]))

    return var_names

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
