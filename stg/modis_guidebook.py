#!/usr/bin/env python
# encoding: utf-8
"""
Provide information about products the system expects in MODIS files, how
they're stored in the files, and how to manage them during the process of
regridding.

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
import os
import logging
import numpy
from datetime import datetime

LOG = logging.getLogger(__name__)

# these are constants for separating the cloud top pressure in to low/mid/high
#       high pressure is < 440
#       low  pressure is >= 680
#       mid pressure is everything in between
HIGH_CLOUD_TOP_PRESSURE_CONST = 440
LOW_CLOUD_TOP_PRESSURE_CONST  = 680

# constants for separating the cloud effective emissivity into thin/thick/opaque
#       thin is < 0.5
#       opaque is >= 0.95
#       thick is everything in between
THIN_CLOUDS_CUTOFF_CONST   = 0.5
OPAQUE_CLOUDS_CUTOFF_CONST = 0.95

# the expected number of files per day if nothing goes wrong
EXPECTED_FILES_PER_DAY      = 288

# variable names expected in the files
CLOUD_PHASE_NAME            = 'Cloud_Phase_Infrared'
CLOUD_TOP_TEMP_NAME         = 'Cloud_Top_Temperature'
CLOUD_TOP_PRESS_NAME        = 'Cloud_Top_Pressure'
CLOUD_TOP_HEIGHT_NAME       = 'Cloud_Top_Height'
CLOUD_EFF_EMISS_NAME        = 'Cloud_Effective_Emissivity'
CLOUD_TOP_PRESS_1KM_NAME    = 'cloud_top_pressure_1km'
CLOUD_TOP_TEMP_1KM_NAME     = 'cloud_top_temperature_1km'
CLOUD_TOP_HEIGHT_1KM_NAME   = 'cloud_top_height_1km'
CLOUD_EFF_EMISS_1KM_NAME    = 'cloud_emissivity_1km'
CLOUD_PHASE_1KM_NAME        = 'Cloud_Phase_Infrared_1km'
CLOUD_EFF_RADIUS_16_NAME    = 'Cloud_Effective_Radius_16'
CLOUD_EFF_RADIUS_37_NAME    = 'Cloud_Effective_Radius_37'
CLOUD_OPTICAL_THICK_NAME    = 'Cloud_Optical_Thickness'
CLOUD_WATER_PATH_NAME       = 'Cloud_Water_Path'
QA_1KM_NAME                 = 'Quality_Assurance_1km'
LATITUDE_NAME               = 'Latitude'
LONGITUDE_NAME              = 'Longitude'
SOLAR_ZENITH_NAME           = 'Solar_Zenith'
SENSOR_ZENITH_NAME          = 'Sensor_Zenith'
SCAN_LINE_TIME_NAME         = 'Scan_Start_Time'
CLOUD_MULTI_LAYER_FLAG_NAME = 'Cloud_Multi_Layer_Flag'
RADIANCE_VARIANCE_NAME      = 'Radiance_Variance'
BRIGHTNESS_TEMP_NAME        = 'Brightness_Temperature'
CLOUD_FRACTION_NAME         = 'Cloud_Fraction'

# TODO, sort out how to differentiate the Cloud Effective Radius
# TODO, currently the fortran is testing:
"""
   if (versionnum<=4 & ~strncmp(filetype,'MAC',3))
      varInFile = 'Effective_Particle_Radius';
   else
      varInFile = 'Cloud_Effective_Radius';
   end
"""
CLOUD_EFF_RADIUS_NAME       = ('Effective_Particle_Radius', 'Cloud_Effective_Radius')

# a list of our navigation variables
NAVIGATION_VAR_NAMES        = [LONGITUDE_NAME,
                               LATITUDE_NAME,
                               SOLAR_ZENITH_NAME,
                               SENSOR_ZENITH_NAME,
                               SCAN_LINE_TIME_NAME,
                               CLOUD_TOP_PRESS_NAME]

# important attribute names
SCALE_ATTR_NAME             = 'scale_factor'
ADD_OFFSET_ATTR_NAME        = 'add_offset'
FILL_VALUE_ATTR_NAME        = '_fillvalue'

# a map of what the caller calls variables to variable names in files
CALLER_VARIABLE_MAP     = {
                            'cloud phase':                  CLOUD_PHASE_NAME,
                            'cloud_phase_infrared' :        CLOUD_PHASE_NAME,
                            'irphase':                      CLOUD_PHASE_NAME,
                            
                            'cloud top temperature':        CLOUD_TOP_TEMP_NAME,
                            'ctt5km':                       CLOUD_TOP_TEMP_NAME,
                            'cloud_top_temperature':        CLOUD_TOP_TEMP_NAME,
                            
                            'cloud top pressure':           CLOUD_TOP_PRESS_NAME,
                            'ctp5km':                       CLOUD_TOP_PRESS_NAME,
                            'pressure':                     CLOUD_TOP_PRESS_NAME,
                            
                            'effective cloud emissivity':   CLOUD_EFF_EMISS_NAME,
                            
                            """ TODO, this won't work for cloud effective radius
                            're':                           CLOUD_EFF_RADIUS_NAME,
                            'effective_particle_radius':    CLOUD_EFF_RADIUS_NAME,
                            'effective radius':             CLOUD_EFF_RADIUS_NAME,
                            'cloud_effective_radius':       CLOUD_EFF_RADIUS_NAME,
                            """
                            
                            'ctp1km':                       CLOUD_TOP_PRESS_1KM_NAME,
                            
                            'ctt1km':                       CLOUD_TOP_TEMP_1KM_NAME,
                            
                            're16':                         CLOUD_EFF_RADIUS_16_NAME,
                            'cloud_effective_radius_16':    CLOUD_EFF_RADIUS_16_NAME,
                            
                            're37':                         CLOUD_EFF_RADIUS_37_NAME,
                            'cloud_effective_radius_37':    CLOUD_EFF_RADIUS_37_NAME,
                            
                            'tau':                          CLOUD_OPTICAL_THICK_NAME,
                            'cloud_optical_thickness':      CLOUD_OPTICAL_THICK_NAME,
                            'optical thickness':            CLOUD_OPTICAL_THICK_NAME,
                            
                            'cwp':                          CLOUD_WATER_PATH_NAME,
                            'cloud_water_path':             CLOUD_WATER_PATH_NAME,
                            'cloud water path':             CLOUD_WATER_PATH_NAME,
                            
                            'retphase':                     QA_1KM_NAME,
                            'retrieval phase':              QA_1KM_NAME,
                            'retqflag':                     QA_1KM_NAME,
                            'retrieval_quality_flag':       QA_1KM_NAME,
                            
                            'latitude':                     LATITUDE_NAME,
                            'lat':                          LATITUDE_NAME,
                            
                            'longitude':                    LONGITUDE_NAME,
                            'lon':                          LONGITUDE_NAME,
                            
                            'solar zenith':                 SOLAR_ZENITH_NAME,
                            'sunzen':                       SOLAR_ZENITH_NAME,
                            
                            'sensor zenith':                SENSOR_ZENITH_NAME,
                            'satzen':                       SENSOR_ZENITH_NAME,
                            'viewing zenith':               SENSOR_ZENITH_NAME,
                            
                            'multilayer':                   CLOUD_MULTI_LAYER_FLAG_NAME,
                            'cloud_multi_layer_flag':       CLOUD_MULTI_LAYER_FLAG_NAME,
                            'overlap':                      CLOUD_MULTI_LAYER_FLAG_NAME,
                            
                            'variance':                     RADIANCE_VARIANCE_NAME,
                            'radiance_variance':            RADIANCE_VARIANCE_NAME,
                            'radiance_variability':         RADIANCE_VARIANCE_NAME,
                            
                            'brightness_temperature':       BRIGHTNESS_TEMP_NAME,
                            'bt':                           BRIGHTNESS_TEMP_NAME,
                          }

DATA_TYPE_TO_USE = { # TODO, eventually differentiate the data type by variable
                    CLOUD_PHASE_NAME:               numpy.float32,
                    CLOUD_TOP_TEMP_NAME:            numpy.float32,
                    CLOUD_TOP_PRESS_NAME:           numpy.float32,
                    CLOUD_TOP_HEIGHT_NAME:          numpy.float32,
                    CLOUD_EFF_EMISS_NAME:           numpy.float32,
                    CLOUD_TOP_PRESS_1KM_NAME:       numpy.float32,
                    CLOUD_TOP_TEMP_1KM_NAME:        numpy.float32,
                    CLOUD_TOP_HEIGHT_1KM_NAME:      numpy.float32,
                    CLOUD_EFF_EMISS_1KM_NAME:       numpy.float32,
                    CLOUD_PHASE_1KM_NAME:           numpy.float32,
                    CLOUD_EFF_RADIUS_16_NAME:       numpy.float32,
                    CLOUD_EFF_RADIUS_37_NAME:       numpy.float32,
                    CLOUD_OPTICAL_THICK_NAME:       numpy.float32,
                    CLOUD_WATER_PATH_NAME:          numpy.float32,
                    QA_1KM_NAME:                    numpy.float32,
                    LATITUDE_NAME:                  numpy.float32,
                    LONGITUDE_NAME:                 numpy.float32,
                    SOLAR_ZENITH_NAME:              numpy.float32,
                    SENSOR_ZENITH_NAME:             numpy.float32,
                    CLOUD_MULTI_LAYER_FLAG_NAME:    numpy.float32,
                    RADIANCE_VARIANCE_NAME:         numpy.float32,
                    BRIGHTNESS_TEMP_NAME:           numpy.float32,
                    CLOUD_FRACTION_NAME:            numpy.float32,
                   }

# a list of the default variables expected in a file, used when no variables are selected by the caller
EXPECTED_VARIABLES_IN_FILE = {
                                CLOUD_TOP_PRESS_NAME,
                                CLOUD_TOP_HEIGHT_NAME,
                                CLOUD_EFF_EMISS_NAME,
                                CLOUD_PHASE_NAME,
                                CLOUD_FRACTION_NAME,
                                # TODO, once the algorithm can handle the 1km data, also add:
                                # CLOUD_TOP_PRESS_1KM_NAME,
                                # CLOUD_TOP_HEIGHT_1KM_NAME,
                                # CLOUD_EFF_EMISS_1KM_NAME,
                                # CLOUD_PHASE_1KM_NAME,

                             }

# TODO, move this up to the general_guidebook
def _clean_off_path_if_needed(file_name_string) :
    """
    remove the path from the file if necessary
    """
    
    return os.path.basename(file_name_string)

def is_MODIS_file (file_name_string) :
    """determine if a file name is the right pattern to represent a MODIS file
    if the file_name_string matches how we expect MODIS files to look return
    TRUE else will return FALSE
    """
    
    temp_name_string = _clean_off_path_if_needed(file_name_string)
    
    return (temp_name_string.startswith('MYD') or temp_name_string.startswith('MOD')) and temp_name_string.endswith('hdf')

def is_MODIS_flat_file (file_name_string) :
    """determine if a flat file is a modis variable
    """

    return ( file_name_string.find(INST_MODIS) >= 0 or
             file_name_string.find(SAT_AQUA)   >= 0 or
             file_name_string.find(SAT_TERRA)  >= 0   )

def parse_datetime_from_filename (file_name_string) :
    """parse the given file_name_string and create an appropriate datetime object
    that represents the datetime indicated by the file name; if the file name does
    not represent a pattern that is understood, None will be returned
    """
    
    temp_name_string = _clean_off_path_if_needed(file_name_string)
    
    datetime_to_return = None
    
    # there are at least two file name formats to parse here
    if temp_name_string.startswith('MYD') or temp_name_string.startswith('MOD') :
        temp = temp_name_string.split('.')
        datetime_to_return = datetime.strptime(temp[1] + temp[2], "A%Y%j%H%M")
        # I confirmed with Nick that this is the correct date format
    
    return datetime_to_return

def get_satellite_from_filename (data_file_name_string) :
    """given a file name, figure out which satellite it's from
    if the file does not represent a known satellite name
    configuration None will be returned
    """
    
    temp_name_string = _clean_off_path_if_needed(data_file_name_string)
    
    satellite_to_return = None
    instrument_to_return = None
    
    if   temp_name_string.startswith("MYD") :
        satellite_to_return  = SAT_AQUA
        instrument_to_return = INST_MODIS
    elif temp_name_string.startswith("MOD") :
        satellite_to_return  = SAT_TERRA
        instrument_to_return = INST_MODIS
    
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

def mask_variable_for_time_gridding (variable_name, variable_data) :
    """given a variable name, return one or more masked areas
    
    the masked areas generated by this function should be treated as
    different variables for the purpose of time gridding
    """
    
    to_return = { }

    # separate the cloud top pressure categories into high/mid/low
    if variable_name == CLOUD_TOP_PRESS_NAME :
        
        to_return[HIGH_MODIFIER] =  variable_data <  HIGH_CLOUD_TOP_PRESSURE_CONST
        to_return[MID_MODIFIER]  = (variable_data >= HIGH_CLOUD_TOP_PRESSURE_CONST) & \
                                   (variable_data <  LOW_CLOUD_TOP_PRESSURE_CONST)
        to_return[LOW_MODIFIER]  =  variable_data >= LOW_CLOUD_TOP_PRESSURE_CONST

    # separate the cloud effective emissivity categories into thin/thick/opaque
    elif variable_name == CLOUD_EFF_EMISS_NAME :

        to_return[THIN_MODIFIER]   =  variable_data <  THIN_CLOUDS_CUTOFF_CONST
        to_return[THICK_MODIFIER]  = (variable_data >= THIN_CLOUDS_CUTOFF_CONST) & \
                                     (variable_data <  OPAQUE_CLOUDS_CUTOFF_CONST)
        to_return[OPAQUE_MODIFIER] =  variable_data >= OPAQUE_CLOUDS_CUTOFF_CONST

    # TODO, separate the cloud phase in to ice/water/unknown categories

    else :
        to_return[""] = numpy.ones(variable_data.shape, dtype=numpy.bool)
    
    return to_return

def get_variable_name_from_flat_file (flat_file_name) :
    """given an fbf file name, figure out which variable is in it
    """
    
    name_to_return = None
    
    for possible_var_name in DATA_TYPE_TO_USE.keys() :
        if flat_file_name.find(possible_var_name) >= 0 :
            name_to_return = possible_var_name
    
    return name_to_return

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
