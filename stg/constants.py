#!/usr/bin/env python
# encoding: utf-8
"""
Module to store all constants.  Any constant needed in more than one
component, or potentially more than one part should be defined here.

Rules/Preferences:
    - All values lowercase
    - strings
    - user-legible (assume that they may be printed in log messages)
    - use == for comparison (not 'is' or 'not' or other)

:author:       Eva Schiffer (evas)
:contact:      eva.schiffer@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
:date:         Jan 2014
:license:      GNU GPLv3

Copyright (C) 2014 - 2015 Space Science and Engineering Center (SSEC),
 University of Wisconsin-Madison.
"""
__docformat__ = "restructuredtext en"

# some constants for time conversion
SECONDS_PER_HOUR          = 3600.0
DEGREES_LON_PER_HOUR      = 15.0
HOURS_PER_DEGREE_LON      = 24.0 / 360.0
HOURS_PER_DAY             = 24.0

# constants for the satellite names
SAT_NOAA_14               = 'noaa-14'
SAT_NOAA_15               = 'noaa-15'
SAT_NOAA_16               = 'noaa-16'
SAT_NOAA_17               = 'noaa-17'
SAT_NOAA_18               = 'noaa-18'
SAT_NOAA_19               = 'noaa-19'
SAT_METOP_A               = 'metop-a'
SAT_METOP_B               = 'metop-b'
SAT_AQUA                  = "aqua"
SAT_TERRA                 = "terra"
ALL_SATS                  = {SAT_NOAA_14,
                             SAT_NOAA_15,
                             SAT_NOAA_16,
                             SAT_NOAA_17,
                             SAT_NOAA_18,
                             SAT_NOAA_19,
                             SAT_METOP_A,
                             SAT_METOP_B,
                             SAT_AQUA,
                             SAT_TERRA}

# constants for the instrument types
INST_AVHRR                = 'avhrr'
INST_HIRS                 = 'hirs'
INST_MODIS                = "modis"

# keys for organizing data
LON_KEY                   = "longitude"
LAT_KEY                   = "latitude"
SCAN_LINE_TIME_KEY        = "scanline-time"
SENSOR_ZENITH_ANGLE_KEY   = "sensor-zenith-angle"
DAY_MASK_KEY              = "day-mask"
NIGHT_MASK_KEY            = "night-mask"
SET_MASK_KEY              = "mask"
# for keeping track of data sets for time gridding
SPACE_GRID_KEY            = "space gridded data"
BLANK_STEM_KEY            = "stem"

# these are categories of separated data used in the system (several may stack)
TEMP_SUFFIX_KEY           = "temp"
DAILY_SPACE_SUFFIX_KEY    = "daily-space"
DAILY_TIME_SUFFIX_KEY     = "daily-time"
MULTI_TIME_SUFFIX_KEY     = "multiday-time"
NOBS_LUT_SUFFIX           = "lut-nobs" # number of observations look up table files
# these represent the categories data may be separated into for space or time gridding
HIGH_MODIFIER             = "high"      # high/mid/low are intended to be used for cloud top pressure categories
MID_MODIFIER              = "mid"
LOW_MODIFIER              = "low"
THIN_MODIFIER             = "thin"      # thin/thick/opaque are intended for use with cloud effective emissivity
THICK_MODIFIER            = "thick"
OPAQUE_MODIFIER           = "opaque"
DAY_SET_KEY               = "day-time"  # day/night/all are intended to separate data from different times of day
NIGHT_SET_KEY             = "night-time"
ALL_SET_KEY               = "all-time"
MORNING_SET_KEY           = "morning-time"
AFTERNOON_SET_KEY         = "afternoon-time"
EVENING_SET_KEY           = "evening-time"
TIME_SETS                 = {DAY_SET_KEY,
                             NIGHT_SET_KEY,
                             ALL_SET_KEY,
                             MORNING_SET_KEY,
                             AFTERNOON_SET_KEY,
                             EVENING_SET_KEY}
# these represent additional data sub-types in both temporary and final files
DENSITY_SUFFIX            = "density"
NOBS_SUFFIX               = "num-observations" # number of observations
NUM_MES_SUFFIX            = "num-measurements" # number of measurements
MEAN_SUFFIX               = "mean"
STD_SUFFIX                = "std"
CLOUD_FRACTION_SUFFIX     = "cloud-fraction" # the cloud fraction
UNCERTAINTY_SUFFIX        = "uncertainty"

# some constants to describe file types
DAILY_SPACE_TYPE          = "daily-space-gridded"
DAILY_TIME_TYPE           = "daily-time-gridded"
MULTIDAY_TIME_TYPE        = "multiday-time-gridded"
NOBS_LUT_TYPE             = "nobs-look-up-table"
ALL_STG_FILE_TYPES        = {DAILY_SPACE_TYPE, DAILY_TIME_TYPE, MULTIDAY_TIME_TYPE, NOBS_LUT_TYPE}

# constants about our output file format
NETCDF_SUFFIX             = "nc"

# the next three are for classifying file sets
TEMP_FILE_TYPE            = "temp-file"
FINAL_OUT_FILE_TYPE       = "final-output-file"
ALL_FILES_TYPE            = "all-files"

# for use in the plotting tool
PLOT_SUFFIX               = "plot_binary"
