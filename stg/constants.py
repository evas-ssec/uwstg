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

INST_AVHRR      = 'avhrr'
INST_HIRS       = 'hirs'
INST_MODIS      = "modis"

LON_KEY         = "longitude"
LAT_KEY         = "latitude"
DAY_MASK_KEY    = "day_mask"
NIGHT_MASK_KEY  = "night_mask"