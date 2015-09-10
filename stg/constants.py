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

Copyright (C) 2014 Space Science and Engineering Center (SSEC),
 University of Wisconsin-Madison.
"""
__docformat__ = "restructuredtext en"

# constants for the satellite names
SAT_NOAA_14     = 'noaa-14'
SAT_NOAA_15     = 'noaa-15'
SAT_NOAA_16     = 'noaa-16'
SAT_NOAA_17     = 'noaa-17'
SAT_NOAA_18     = 'noaa-18'
SAT_NOAA_19     = 'noaa-19'
SAT_METOP_A     = 'metop-a'
SAT_METOP_B     = 'metop-b'
SAT_AQUA        = "aqua"
SAT_TERRA       = "terra"
ALL_SATS        = {SAT_NOAA_14,
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
INST_AVHRR      = 'avhrr'
INST_HIRS       = 'hirs'
INST_MODIS      = "modis"

LON_KEY         = "longitude"
LAT_KEY         = "latitude"
DAY_MASK_KEY    = "day_mask"
NIGHT_MASK_KEY  = "night_mask"

PLOT_SUFFIX          = "plot_binary"

# these are suffixes used for temporary files
DAY_TEMP_SUFFIX           = "daytemp"
NIGHT_TEMP_SUFFIX         = "nighttemp"
DAY_DENSITY_TEMP_SUFFIX   = "daydensitytemp"
NIGHT_DENSITY_TEMP_SUFFIX = "nightdensitytemp"
DAY_NOBS_TEMP_SUFFIX      = "daynobstemp"
NIGHT_NOBS_TEMP_SUFFIX    = "nightnobstemp"
ALL_TEMP_SUFFIX           = "alltemp"
ALL_DENSITY_TEMP_SUFFIX   = "alldensitytemp"
ALL_NOBS_TEMP_SUFFIX      = "allnobstemp"

# these are suffixes used for the final, packed daily space gridded files
DAY_SUFFIX                = "dayfinal"
NIGHT_SUFFIX              = "nightfinal"
DAY_NOBS_SUFFIX           = "daynobsfinal"
NIGHT_NOBS_SUFFIX         = "nightnobsfinal"
DAY_NUM_MES_SUFFIX        = "daynummeasurements"
NIGHT_NUM_MES_SUFFIX      = "nightnummeasurements"
ALL_SUFFIX                = "allfinal"
ALL_NOBS_SUFFIX           = "allnobsfinal"
ALL_NUM_MES_SUFFIX        = "allnummeasurements"
DAILY_NOBS_KEY            = "nobsfinal"
DAILY_SPACE_SUFFIXES      = {DAY_SUFFIX,
                             NIGHT_SUFFIX,
                             DAY_NOBS_SUFFIX,
                             NIGHT_NOBS_SUFFIX,
                             DAY_NUM_MES_SUFFIX,
                             NIGHT_NUM_MES_SUFFIX,
                             ALL_SUFFIX,
                             ALL_NOBS_SUFFIX,
                             ALL_NUM_MES_SUFFIX,}

# these represent the categories data may be separated into for time gridding
HIGH_MODIFIER             = "high"      # high/mid/low are intended to be used for cloud top pressure categories
MID_MODIFIER              = "mid"
LOW_MODIFIER              = "low"

# these are suffixes for time gridded files
DAILY_MEAN_SUFFIX         = "daily-mean"
DAILY_NUM_MES_SUFFIX      = "daily-num-measurements"
DAILY_NOBS_SUFFIX         = "daily-num-observations"
DAILY_MIN_SUFFIX          = "daily-min"
DAILY_MAX_SUFFIX          = "daily-max"
DAILY_STD_SUFFIX          = "daily-std"
DAILY_FRACTION_SUFFIX     = "daily-cloud-fraction"
DAILY_UNCERTAINTY_SUFFIX  = "daily-uncertainty"
DAILY_TIME_SUFFIXES       = {DAILY_MEAN_SUFFIX,
                             DAILY_NUM_MES_SUFFIX,
                             DAILY_NOBS_SUFFIX,
                             DAILY_MIN_SUFFIX,
                             DAILY_MAX_SUFFIX,
                             DAILY_STD_SUFFIX,
                             DAILY_FRACTION_SUFFIX,
                             DAILY_UNCERTAINTY_SUFFIX,}

# look up table files
NOBS_LUT_SUFFIX           = "lut-nobs"

# some constants to describe file types
DAILY_SPACE_TYPE          = "daily-space-gridded-file"
DAILY_TIME_TYPE           = "daily-time-gridded-file"
MULTIDAY_TIME_TYPE        = "multi-day-time-gridded-file"
NOBS_LUT_TYPE             = "num-obs-look-up-table-file"

# keys for keeping track of data sets
DAY_SET_KEY               = "day"
NIGHT_SET_KEY             = "night"
ALL_SET_KEY               = "all"
DAILY_SET                 = "daily"
ALL_SET                   = "all"
SET_MASK_KEY              = "mask"
SET_TEMP_DENSITY_SUFF_KEY = "temp density suffix"
SET_TEMP_NOBS_SUFF_KEY    = "temp nobs suffix"
SET_TEMP_DATA_SUFF_KEY    = "temp data suffix"
SET_FINAL_DATA_SUFF_KEY   = "output data suffix"
SET_FINAL_NOBS_SUFF_KEY   = "output nobs suffix"
SET_FINAL_NMES_SUFF_KEY   = "output num measurements suffix"

# for keeping track of data sets for time gridding
NOBS_KEY                  = "nobs"
SPACE_GRID_KEY            = "space gridded data"
BLANK_STEM_KEY            = "stem"
