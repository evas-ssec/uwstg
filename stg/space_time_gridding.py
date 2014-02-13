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
        desired_variables = list(args) if len(args) > 0 else [ ]
        input_path        = options.inputPath
        output_path       = options.outputPath
        min_scan_angle    = options.minScanAngle
        grid_degrees      = float(options.gridDegrees)
        
        # determine the grid size in number of elements
        grid_lon_size    = int(math.ceil(360.0 / grid_degrees))
        grid_lat_size    = int(math.ceil(180.0 / grid_degrees))
        space_grid_shape = (grid_lon_size, grid_lat_size) # TODO, is this the correct order?
        
        # look through our files and figure out what variables we expect from them
        possible_files    = os.listdir(input_path)
        expected_vars     = { }
        all_vars          = set()
        date_time_temp    = None
        for file_name in sorted(possible_files) :
            expected_vars[file_name] = general_guidebook.get_variable_names (file_name, user_requested_names=desired_variables)
            # if this file has no variables, remove it from our files for consideration
            if len(expected_vars[file_name]) <= 0 :
                del expected_vars[file_name]
                possible_files.remove(file_name)
            # otherwise, add the variables we found to our list of all variables and try to get a time from the file
            else :
                all_vars.update(expected_vars[file_name])
                date_time_temp = general_guidebook.parse_datetime_from_filename(file_name) if date_time_temp is None else date_time_temp
        
        # check to make sure our intermediate file names don't exist already
        for var_name in all_vars :
            
            for suffix in io_manager.ALL_EXPECTED_SUFFIXES :
                # TODO, pull satellite and algorithm too
                temp_stem = io_manager.build_name_stem(var_name, date_time=date_time_temp, satellite=None, algorithm=None, suffix=suffix)
                temp_name = fbf.filename(temp_stem, TEMP_DATA_TYPE, shape=(space_grid_shape))
                if os.path.exists(os.path.join(output_path, temp_name)) :
                    LOG.warn ("Cannot process files because matching temporary or output files exist in the output directory.")
                    return
        
        # loop to deal with data from each of the files
        failed_files    = 0
        sucessful_files = 0
        for each_file in sorted(possible_files) :
            
            full_file_path = os.path.join(input_path, each_file)
            
            LOG.debug("Processing file: " + full_file_path)
            
            # load the aux data
            file_object, temp_aux_data = io_manager.load_aux_data(full_file_path,
                                                                  min_scan_angle)
            
            ok_file = True
            try :
                # calculate the indecies for the space grid based on the aux data
                # (we can do this now since the lon/lat is the same for each variable in the file)
                day_lon_index, day_lat_index, night_lon_index, night_lat_index = space_gridding.calculate_index_from_nav_data(temp_aux_data,
                                                                                                                              grid_degrees)
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
                    
                    # split the variable by day/night
                    day_var_data   = var_data[temp_aux_data[DAY_MASK_KEY]]
                    night_var_data = var_data[temp_aux_data[NIGHT_MASK_KEY]]
                    
                    
                    ok_file = True
                    try :
                        # space grid the data using the indexes we calculated earlier
                        day_space_grid,   day_density_map,   day_nobs,   day_max_depth   = space_gridding.space_grid_data (grid_lon_size, grid_lat_size,
                                                                                                                           day_var_data,
                                                                                                                           day_lon_index, day_lat_index)
                        night_space_grid, night_density_map, night_nobs, night_max_depth = space_gridding.space_grid_data (grid_lon_size, grid_lat_size,
                                                                                                                           night_var_data,
                                                                                                                           night_lon_index, night_lat_index)
                    except Exception, e :
                    
                        LOG.warn("Unable to process variable data space gridding for file: " + full_file_path)
                        LOG.warn("This variable will not be processed.")
                        
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        LOG.debug(traceback.format_exception(exc_type, exc_value, exc_traceback))
                        
                        ok_file       = False
                        failed_files += 1
                    
                    if ok_file :
                        
                        # save the space grids and density info for this variable and it's density map to files
                        
                        # day related files
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.DAY_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, day_space_grid, TEMP_DATA_TYPE)
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.DAY_DENSITY_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, day_density_map, TEMP_DATA_TYPE)
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.DAY_NOBS_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, day_nobs, TEMP_DATA_TYPE)
                        
                        # night related files
                        
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.NIGHT_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, night_space_grid, TEMP_DATA_TYPE)
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.NIGHT_DENSITY_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, night_density_map, TEMP_DATA_TYPE)
                        io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                 satellite=None, algorithm=None,
                                                                                 suffix=io_manager.NIGHT_NOBS_TEMP_SUFFIX),
                                                     space_grid_shape, output_path, night_nobs, TEMP_DATA_TYPE)
                
            # make sure each file is closed when we're done with it
            io_manager.close_file(full_file_path, file_object)
            
            # if we got to here we processed the file correctly
            sucessful_files += 1
        
        # collapse the per variable space grids to remove excess NaNs
        for variable_name in all_vars :
            
            LOG.debug("Packing space data for variable: " + variable_name)
            
            # load the variable's density maps
            var_workspace     = Workspace.Workspace(dir=output_path)
            day_var_density   = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                         satellite=None, algorithm=None,
                                                                         suffix=io_manager.DAY_DENSITY_TEMP_SUFFIX)][:] 
            night_var_density = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                         satellite=None, algorithm=None,
                                                                         suffix=io_manager.NIGHT_DENSITY_TEMP_SUFFIX)][:] 
            
            # only do the day data if we have some
            if numpy.sum(day_var_density) > 0 :
                
                # load the sparse space grid
                day_var_data      = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                             satellite=None, algorithm=None,
                                                                             suffix=io_manager.DAY_TEMP_SUFFIX)][:]
                
                # collapse the space grid
                final_day_data    = space_gridding.pack_space_grid(day_var_data,   day_var_density)
                
                # save the final array to an appropriately named file
                io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                        satellite=None, algorithm=None,
                                                                        suffix=io_manager.DAY_SUFFIX),
                                             space_grid_shape, output_path, final_day_data,
                                             TEMP_DATA_TYPE, file_permissions="w")
                
                # load the nobs file
                nobs_counts       = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                             satellite=None, algorithm=None,
                                                                             suffix=io_manager.DAY_NOBS_TEMP_SUFFIX)][:]
                
                # collapse the nobs
                nobs_final        = numpy.sum(nobs_counts, axis=0)
                
                # save the final nobs array to an appropriately named file
                io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                        satellite=None, algorithm=None,
                                                                        suffix=io_manager.DAY_NOBS_SUFFIX),
                                             space_grid_shape, output_path,
                                             nobs_final, TEMP_DATA_TYPE, file_permissions="w")
                
            else :
                LOG.warn("No day data was found for variable " + variable_name + ". Day files will not be written.")
            
            # only do night data if we have some
            if numpy.sum(night_var_density) > 0 :
                
                # load the sparse space grid
                night_var_data      = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                             satellite=None, algorithm=None,
                                                                             suffix=io_manager.NIGHT_TEMP_SUFFIX)][:]
                
                # collapse the space grid
                final_night_data    = space_gridding.pack_space_grid(night_var_data,   night_var_density)
                
                # save the final array to an appropriately named file
                io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                        satellite=None, algorithm=None,
                                                                        suffix=io_manager.NIGHT_SUFFIX),
                                             space_grid_shape, output_path, final_night_data,
                                             TEMP_DATA_TYPE, file_permissions="w")
                
                # load the nobs file
                nobs_counts       = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                             satellite=None, algorithm=None,
                                                                             suffix=io_manager.NIGHT_NOBS_TEMP_SUFFIX)][:]
                
                # collapse the nobs
                nobs_final        = numpy.sum(nobs_counts, axis=0)
                
                # save the final nobs array to an appropriately named file
                io_manager.save_data_to_file(io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                        satellite=None, algorithm=None,
                                                                        suffix=io_manager.NIGHT_NOBS_SUFFIX),
                                             space_grid_shape, output_path,
                                             nobs_final, TEMP_DATA_TYPE, file_permissions="w")
                
            else :
                LOG.warn("No night data was found for variable " + variable_name + ". Night files will not be written.")
        
        LOG.debug("Successfully processed " + str(sucessful_files) + " files and failed to process " + str(failed_files) + " files for this day.")
        
        # remove the extra temporary files in the output directory
        remove_suffixes = ["*" + p + "*" for p in io_manager.EXPECTED_TEMP_SUFFIXES]
        remove_file_patterns(output_path, remove_suffixes)
    
    def stats_day(*args) :
        """given files of daily space gridded data, calculate daily stats
        given an input directory that contains appropriate files,
        calculate daily stats and put the resulting gridded files
        for that day in the output directory.
        
        Note: the output directory will also be used for intermediary working
        files.
        """
        
        # set up some of our input from the caller for easy access
        desired_variables = list(args) if len(args) > 0 else [ ]
        input_path        = options.inputPath
        output_path       = options.outputPath
        min_scan_angle    = options.minScanAngle
        grid_degrees      = float(options.gridDegrees)
    
    def stats_month(*args) :
        """given a month of daily space gridded data, calculate montly stats
        given an input directory that contains appropriate files,
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

