#!/usr/bin/env python
# encoding: utf-8
"""
This module handles statistics across multiple times of pre-space-gridded data.

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

import numpy

def create_sample_size_cutoff_mask (data_array, nobs_array,
                                    overall_nobs_array,
                                    fixed_cutoff=None, dynamic_std_cutoff=None) :
    """
    given a data array, a matching number of observations array, and an overall number of observations
    for the full time period, create a mask of which data cells should be discarded from this data set.
    
    Cutoffs can either be a fixed number (ie. fixed_cutoff=2 would cause any cell with two or less
    observations to be excluded) or a dynamic standard deviation range (ie. dynamic_std_cutoff=1.0
    would cause any observations more than one standard deviation less than the overall mean for that
    cell to be excluded).
    """
    
    # start out by assuming all our data is good
    bad_data_mask = numpy.zeros(nobs_array.shape, dtype=numpy.bool)
    
    # if we have a fixed cutoff, apply that
    if fixed_cutoff is not None :
        bad_data_mask[nobs_array <= fixed_cutoff] = True
    
    # if we have a dynamic cutoff, apply that
    if dynamic_std_cutoff is not None :
        
        mean_overall_nobs = numpy.mean(overall_nobs_array, axis=0)
        std_overall_nobs  = numpy.std (overall_nobs_array, axis=0)
        cutoff_values     = mean_overall_nobs - (std_overall_nobs * dynamic_std_cutoff)
        
        bad_data_mask[nobs_array < cutoff_values] = True
    
    return bad_data_mask

def calculate_partial_weighted_time_average (daily_data_array, daily_nobs_array,
                                             overall_sum_num_measurments_array, overall_sum_nobs_array) :
    """
    given a set of daily data and their matching number of observations array as well as the
    sum of the overall number of measurments and number of observations per cell (ie. the collapsed
    2D arrays summed over the whole time period being analyzed), calculate the contribution from this
    day to the weighted time average
    
    Note: to calculate the weighted time average over the full period, add up all of the daily
    partial weighted time averages in that period
    """
    
    good_measurments_this_day = numpy.sum(numpy.isfinite(daily_data_array), axis=0)
    this_day_weighted_average = ( (daily_data_array * (good_measurments_this_day / daily_nobs_array))
                                / (overall_sum_num_measurments_array / overall_sum_nobs_array) )
    
    return this_day_weighted_average






