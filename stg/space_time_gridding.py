#!/usr/bin/env python
# encoding: utf-8
"""
The main module for the Space Time Gridding application. This module handles
the command line input from the user and coordinates the other modules.

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

import logging
import pkg_resources
import os
import math
import glob
import traceback
import sys

import numpy

from collections import defaultdict

import keoni.fbf.workspace as Workspace
import keoni.fbf       as fbf

from stg.constants import *
import stg.general_guidebook as general_guidebook
import stg.io_manager        as io_manager
import stg.space_gridding    as space_gridding

# TODO, in the long run handle the dtype more flexibly
TEMP_DATA_TYPE = numpy.dtype(numpy.float32)

# TODO, need to confirm with Nadia that this is the cutoff she wants
EXPECTED_FRACTION_OF_FILES_PER_DAY = 2.0 / 3.0

LOG = logging.getLogger(__name__)

def get_version_string() :
    version_num = pkg_resources.require('spacetimegrid')[0].version
    
    return "Space Time Gridding, version " + str(version_num) 

def remove_file_patterns(path, *args):
    """Remove files that were created from a previous run,
    including temporary files, that may conflict with
    future processing.

    """
    for pat_list in args:
        for pat in pat_list:
            full_pat = os.path.join(path, pat)
            for f in glob.glob(full_pat):
                _safe_remove(f)

def _safe_remove(fn):
    """Remove the file `fn` if you can, if not log an error message,
    but continue on.

    :Parameters:
        fn : str
            Filename of the file to be removed.
    """
    try:
        LOG.info("Removing %s" % fn)
        os.remove(fn)
    except StandardError:
        LOG.error("Could not remove %s" % fn)

def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
examples:

python -m space_time_gridding 

"""
    
    # set the options available to the user on the command line
    parser = optparse.OptionParser(usage)
    
    # message related options
    parser.add_option('-q', '--quiet', dest="quiet",
                    action="store_true", default=False, help="only error output")
    parser.add_option('-v', '--verbose', dest="verbose",
                    action="store_true", default=False, help="enable more informational output")   
    parser.add_option('-w', '--debug', dest="debug",
                    action="store_true", default=False, help="enable debug output")
    parser.add_option('-n', '--version', dest='version',
                      action="store_true", default=False, help="view the STG version")
    
    # output generation related options
    parser.add_option('-i', '--input',  dest="inputPath",  type='string', default='./',
                      help="set path for the input directory")
    parser.add_option('-o', '--output', dest="outputPath", type='string', default='./out/',
                      help="set path for the output directory")
    
    # options related to space or time gridding
    parser.add_option('-g', '--grid_degrees', dest="gridDegrees", type='float', default=1.0,
                      help="set the size of the output grid's cells in degrees")
    parser.add_option('-a', '--min_scan_angle', dest="minScanAngle", type='float', default=60.0,
                      help="the minimum scan angle that will be considered useful")
    parser.add_option('-d', '--do_process_with_little_data', dest="overrideMinCheck",
                      action="store_true", default=False, help="run the full daily compilation even if many files are missing or unreadable")
    
    # parse the uers options from the command line
    options, args = parser.parse_args()
    
    # set up the logging level based on the options the user selected on the command line
    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)
    
    # display the version
    if options.version :
        print (get_version_string() + '\n')

    commands = {}
    prior = None
    prior = dict(locals())
    
    """
    The following functions represent available menu selections
    """
    
    def space_day(*args) :
        """grid one day of input files in space
        given an input directory that contains appropriate files,
        grid them in space and put the resulting gridded files
        for that day in the output directory.
        
        Note: the output directory will also be used for intermediary working
        files.
        """
        
        # set up some of our input from the caller for easy access
        desired_variables  = list(args) if len(args) > 0 else [ ]
        input_path         = options.inputPath
        output_path        = options.outputPath
        min_scan_angle     = options.minScanAngle
        grid_degrees       = float(options.gridDegrees)
        
        # determine the grid size in number of elements
        grid_lon_size      = int(math.ceil(360.0 / grid_degrees))
        grid_lat_size      = int(math.ceil(180.0 / grid_degrees))
        space_grid_shape   = (grid_lon_size, grid_lat_size) # TODO, is this the correct order?
        
        # look through our files and figure out what variables we expect from them
        possible_files     = os.listdir(input_path)
        expected_vars      = { }
        all_vars           = set()
        date_time_temp     = None
        expected_num_files = None
        satellite          = None
        instrument         = None
        for file_name in sorted(possible_files) :
            expected_vars[file_name] = general_guidebook.get_variable_names (file_name, user_requested_names=desired_variables)
            # if this file has no variables, remove it from our files for consideration
            if len(expected_vars[file_name]) <= 0 :
                del expected_vars[file_name]
                possible_files.remove(file_name)
            # otherwise, add the variables we found to our list of all variables and try to get a time from the file
            else :
                all_vars.update(expected_vars[file_name])
                # if we don't have it yet, update some general information about this run based on the file name
                temp_sat, temp_inst   = general_guidebook.get_satellite_from_filename(file_name)
                satellite             = temp_sat  if satellite  is None else satellite
                instrument            = temp_inst if instrument is None else instrument
                date_time_temp        = general_guidebook.parse_datetime_from_filename(file_name) if date_time_temp     is None else date_time_temp
                expected_num_files    = general_guidebook.get_expected_files_per_day(instrument)  if expected_num_files is None else expected_num_files
        
        # check to make sure our intermediate file names don't exist already
        for var_name in all_vars :
            
            for suffix in io_manager.ALL_EXPECTED_SUFFIXES :
                # TODO, pull algorithm too
                temp_stem = io_manager.build_name_stem(var_name, date_time=date_time_temp, satellite=satellite, algorithm=None, suffix=suffix)
                temp_name = fbf.filename(temp_stem, TEMP_DATA_TYPE, shape=(space_grid_shape))
                if os.path.exists(os.path.join(output_path, temp_name)) :
                    LOG.warn ("Cannot process files because matching temporary or output files exist in the output directory.")
                    return
        
        # loop to deal with data from each of the files
        failed_files       = 0
        sucessful_files    = 0
        abstract_data_sets = io_manager.get_expected_abstract_sets(instrument)
        for each_file in sorted(possible_files) :
            
            full_file_path = os.path.join(input_path, each_file)
            
            LOG.debug("Processing file: " + full_file_path)
            
            # load the aux data
            file_object, temp_aux_data = io_manager.load_aux_data(full_file_path, min_scan_angle)
            # figure out what data sets we need to process
            data_sets = io_manager.get_expected_data_sets_from_aux_data (instrument, temp_aux_data)
            
            ok_file     = True
            lon_indices = { }
            lat_indices = { }
            try :
                
                # calculate the indecies for the space grid based on the navigation data
                # (we can do this now since the lon/lat is the same for each variable in the file)
                for set_key in data_sets.keys() :
                    
                    set_mask      = data_sets[set_key][SET_MASK_KEY]
                    temp_lon_data = data_sets[set_key][LON_KEY][set_mask]
                    temp_lat_data = data_sets[set_key][LAT_KEY][set_mask]
                    
                    lat_index, lon_index = space_gridding.calculate_index_from_nav_data(temp_lat_data, temp_lon_data, grid_degrees)
                    lat_indices[set_key] = lat_index
                    lon_indices[set_key] = lon_index
                
            except Exception, e :
                
                LOG.warn("Unable to process basic space gridding for file: " + full_file_path)
                LOG.warn("This file will not be processed.")
                
                exc_type, exc_value, exc_traceback = sys.exc_info()
                LOG.debug(traceback.format_exception(exc_type, exc_value, exc_traceback))
                
                ok_file       = False
                failed_files += 1
            
            if ok_file :
                
                # loop to load each variable in the file and process it
                for variable_name in expected_vars[each_file] :
                    
                    LOG.debug("Processing variable: " + variable_name)
                    
                    # load the variable
                    file_object, var_data = io_manager.load_variable_from_file (variable_name,
                                                                                file_path=full_file_path,
                                                                                file_object=file_object)
                    
                    # split the variable data by sets
                    separated_data = { }
                    for set_key in data_sets.keys() :
                        
                        separated_data[set_key] = var_data[data_sets[set_key][SET_MASK_KEY]]
                    
                    ok_file = True
                    space_grids  = { }
                    density_maps = { }
                    nobs         = { }
                    max_depths   = { }
                    try :
                        
                        # space grid the data using the indexes we calculated earlier
                        for set_key in data_sets.keys() :
                            temp_space_grid, temp_density_map, temp_nobs, temp_max_depth = space_gridding.space_grid_data(grid_lat_size,
                                                                                                                          grid_lon_size,
                                                                                                                          separated_data[set_key],
                                                                                                                          lat_indices[set_key],
                                                                                                                          lon_indices[set_key])
                            space_grids[set_key]  = temp_space_grid
                            density_maps[set_key] = temp_density_map
                            nobs[set_key]         = temp_nobs
                            max_depths[set_key]   = temp_max_depth
                        
                    except Exception, e :
                    
                        LOG.warn("Unable to process variable data space gridding for file: " + full_file_path)
                        LOG.warn("This variable will not be processed.")
                        
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        LOG.debug(traceback.format_exception(exc_type, exc_value, exc_traceback))
                        
                        ok_file       = False
                        failed_files += 1
                    
                    if ok_file :
                        
                        # save the space grids and density info for this variable and it's density map to files
                        for set_key in data_sets.keys() :
                            
                            # save the gridded data
                            io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                    satellite=satellite, algorithm=None,
                                                                                    suffix=data_sets[set_key][SET_TEMP_DATA_SUFF_KEY]),
                                                        space_grid_shape, output_path, space_grids[set_key], TEMP_DATA_TYPE)
                            # save the grid density map
                            io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                    satellite=satellite, algorithm=None,
                                                                                    suffix=data_sets[set_key][SET_TEMP_DENSITY_SUFF_KEY]),
                                                        space_grid_shape, output_path, density_maps[set_key], TEMP_DATA_TYPE)
                            # save the number of observations grid
                            io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                    satellite=satellite, algorithm=None,
                                                                                    suffix=data_sets[set_key][SET_TEMP_NOBS_SUFF_KEY]),
                                                        space_grid_shape, output_path, nobs[set_key], TEMP_DATA_TYPE)
            
            # make sure each file is closed when we're done with it
            io_manager.close_file(full_file_path, file_object)
            
            # if we got to here we processed the file correctly
            sucessful_files += 1
        
        LOG.debug("Successfully processed " + str(sucessful_files) + " files and failed to process " + str(failed_files) + " files for this day.")
        
        # warn the user if we have fewer files than we need for this instrument
        if sucessful_files < (expected_num_files * EXPECTED_FRACTION_OF_FILES_PER_DAY) :
            LOG.warn("Processed " + str(sucessful_files)    + " files successfully for this day.")
            LOG.warn("Expected  " + str(expected_num_files) + " files for this instrument type.")
            
            if options.overrideMinCheck :
                LOG.warn ("Daily file(s) will be produced, but data may be unusable for this day.")
            else :
                LOG.critical("Daily file(s) will not be produced for this day due to lack of data.")
                LOG.critical("If you wish to produce the daily file(s), rerun the program using the \'-d\' option.")
        
        # only collect the daily data if we have enough files or have turned off the minimum check
        if ( (sucessful_files >= (expected_num_files * EXPECTED_FRACTION_OF_FILES_PER_DAY)) or
             (options.overrideMinCheck) ):
            
            # collapse the per variable space grids to remove excess NaNs
            for variable_name in all_vars :
                
                LOG.debug("Packing space data for variable: " + variable_name)
                
                # load the variable's density maps
                var_workspace     = Workspace.Workspace(dir=output_path)
                
                for set_key in abstract_data_sets.keys()  :
                    
                    # load the density
                    temp_density = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                             satellite=satellite, algorithm=None,
                                                                             suffix=abstract_data_sets[set_key][SET_TEMP_DENSITY_SUFF_KEY])][:]
                    
                    # only process the final data if it exists
                    if numpy.sum(temp_density) > 0 :
                        
                        # load the sparse space grid
                        var_data = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                            satellite=satellite, algorithm=None,
                                                                            suffix=abstract_data_sets[set_key][SET_TEMP_DATA_SUFF_KEY])][:]
                        
                        # collapse the space grid
                        final_data = space_gridding.pack_space_grid(var_data, temp_density)
                        
                        # save the final array to an appropriately named file
                        io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                satellite=satellite, algorithm=None,
                                                                                suffix=abstract_data_sets[set_key][SET_FINAL_DATA_SUFF_KEY]),
                                                     space_grid_shape, output_path, final_data,
                                                     TEMP_DATA_TYPE, file_permissions="w")
                        
                        # load the nobs file
                        nobs_counts = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                     satellite=satellite, algorithm=None,
                                                                                     suffix=abstract_data_sets[set_key][SET_TEMP_NOBS_SUFF_KEY])][:]
                        
                        # collapse the nobs for the whole day
                        nobs_final = numpy.sum(nobs_counts, axis=0)
                        
                        # save the final nobs array to an appropriately named file
                        io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                satellite=satellite, algorithm=None,
                                                                                suffix=abstract_data_sets[set_key][SET_FINAL_NOBS_SUFF_KEY]),
                                                     space_grid_shape, output_path,
                                                     nobs_final, TEMP_DATA_TYPE, file_permissions="w")
                        
                    else :
                        LOG.warn("No " + set_key + " data was found for variable " + variable_name + ". Corresponding files will not be written.")
        
        # remove the extra temporary files in the output directory
        remove_suffixes = ["*" + p + "*" for p in io_manager.EXPECTED_TEMP_SUFFIXES]
        remove_file_patterns(output_path, remove_suffixes)
    
    def stats_day(*args) :
        """given a day worth of files of daily space gridded data, calculate daily stats
        
        given an input directory that contains space gridded files for a day,
        calculate daily stats and put the resulting time gridded files
        for that day in the output directory.
        
        Note: the output directory will also be used for intermediary working
        files.
        """
        
        # set up some of our input from the caller for easy access
        desired_variables = list(args) if len(args) > 0 else [ ]
        input_path        = options.inputPath
        output_path       = options.outputPath
        
        
        # check the directory for sets of daily files
        
        # for each set of daily files
        
        #   load the main data
        #   load the nobs
        
        #   calculate the std, min, max, and (weighted or non-weighted) average
        
        #   save the various stats to files
        
    
    def stats_month(*args) :
        """given a month of daily space gridded data, calculate montly stats
        given an input directory that contains appropriate daily stats,
        calculate monthly stats and put the resulting gridded files
        for that month in the output directory.
        
        Note: the output directory will also be used for intermediary working
        files.
        """
        
        # set up some of our input from the caller for easy access
        desired_variables = list(args) if len(args) > 0 else [ ]
        input_path        = options.inputPath
        output_path       = options.outputPath
        min_scan_angle    = options.minScanAngle
        grid_degrees      = float(options.gridDegrees)
        
    
    
    # all the local public functions are considered part of the application, collect them up
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    # if what the user asked for is not one of our existing functions, print the help
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        help()
        return 9
    else:
        # call the function the user named, given the arguments from the command line  
        rc = locals()[args[0]](*args[1:])
        return 0 if rc is None else rc
    
    return 0 # it shouldn't be possible to get here any longer

if __name__=='__main__':
    sys.exit(main())

