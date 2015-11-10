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

Copyright (C) 2014 - 2015 Space Science and Engineering Center (SSEC),
 University of Wisconsin-Madison.
"""
__docformat__ = "restructuredtext en"

import logging, pkg_resources, os, math, glob, traceback, sys, time
from collections import defaultdict

import numpy

from netCDF4 import Dataset

import keoni.fbf.workspace as Workspace
import keoni.fbf       as fbf

from   stg.constants import *
import stg.general_guidebook as general_guidebook
import stg.io_manager        as io_manager
import stg.space_gridding    as space_gridding
import stg.time_gridding     as time_gridding
import stg.stg_util          as stg_util

# TODO, in the long run handle the dtype more flexibly
TEMP_DATA_TYPE = numpy.dtype(numpy.float32)

# how much extra we will expand arrays by
ARRAY_GROWTH_FACTOR = 1.25

# have confirmed with Nadia that this is the cutoff she wants for now
EXPECTED_FRACTION_OF_FILES_PER_DAY = 2.0 / 3.0

LOG = logging.getLogger(__name__)

def _get_version_string() :
    version_num = pkg_resources.require('spacetimegrid')[0].version
    
    return "Space Time Griding, version " + str(version_num)

def _remove_file_patterns(path, *args):
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

def _expand_array_if_needed (in_array, min_size, fill_value=numpy.nan) :
    """expand the first dimension of the in_array to be at least min_size size, keeping the existing array data

    Note: this currently expects three dimensional numpy arrays
    """

    to_return = in_array

    if in_array.shape[0] < min_size :

        new_shape = (min_size, in_array.shape[1], in_array.shape[2])

        to_return = numpy.ones(new_shape, dtype=in_array.dtype) * fill_value
        to_return[0:in_array.shape[0], :, :] = in_array

    return to_return

def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
examples:

stg space_griding_day -i /input/path -g 0.5 -t
stg time_griding_day -i /input/path -o /output/path
stg time_griding_multiday
stg flatfile_to_netCDF -o /output/path
stg make_nobs_lut -i /input/path

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
                      action="store_true", default=False, help="print the STG version")
    
    # output generation related options
    parser.add_option('-i', '--input',  dest="inputPath",  type='string', default='./',
                      help="set path for the input directory; defaults to ./")
    parser.add_option('-o', '--output', dest="outputPath", type='string', default='./out/',
                      help="set path for the output directory; defaults to ./out/")
    
    # options related to space or time gridding
    parser.add_option('-g', '--grid_degrees', dest="gridDegrees", type='float', default=0.5,
                      help="set the size of the output grid's cells in degrees; defaults to 0.5")
    parser.add_option('-a', '--min_scan_angle', dest="minScanAngle", type='float', default=32.0,
                      help="the minimum scan angle that will be considered useful; defaults to 32.0")
    parser.add_option('-p', '--do_process_with_little_data', dest="overrideMinCheck",
                      action="store_true", default=False, help="run the full daily compilation even if many files are missing or unreadable")
    parser.add_option('-f', '--fixed_nobs_cutoff', dest="fixedNobsCutoff", type='float', default=None,
                      help="the minimum number of nobs that must be present for data to be considered when time gridding; "
                           "by default there will be no check")
    parser.add_option('-d', '--dynamic_nobs_cutoff', dest="dynamicNobsCutoff", type='float', default=None,
                      help="the minimum nobs that must be present for data to be considered when time gridding," +
                           " expressed a fraction of the std from the mean; by default there will be no check")
    parser.add_option('-l', '--nobs_lut', dest="nobsLUT", type='string', default=None,
                      help="a long term look up table for use with the --dynamic_nobs_cutoff option; by default no lut will be loaded")
    parser.add_option('-t', '--day_night_together', dest='keep_day_night_together',
                      action="store_true", default=False, help="instead of separating day and night data, process them together")
    parser.add_option('-m', '--multiple_overpasses_per_cell', dest="allow_multiple_overpasses_per_cell",
                      action="store_true", default=False, help="allow multiple overpasses to fall in the same grid cell; "
                                                               "the default is to only allow the best overpass in each cell"
                                                               " (best is determined by sensor angle)")
    
    # parse the users options from the command line
    options, args = parser.parse_args()
    
    # set up the logging level based on the options the user selected on the command line
    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)
    
    # display the version
    if options.version :
        print (_get_version_string() + '\n')

    commands = {}
    prior = None
    prior = dict(locals())
    
    ##### The following functions represent available menu selections #####
    
    def space_griding_day(*args) :
        """grid one day of input files in space
        given an input directory that contains appropriate files,
        grid them in space and put the resulting gridded files
        for that day in the output directory.
        
        Note: the output directory will also be used for intermediary working
        files.
        """
        
        # set up some of our input from the caller for easy access
        desired_variables  = list(args) if len(args) > 0 else [ ]
        input_path         = stg_util.clean_path(options.inputPath)
        output_path        = stg_util.clean_path(options.outputPath)
        stg_util.setup_dir_if_needed(output_path, "output")
        min_scan_angle     = options.minScanAngle
        grid_degrees       = float(options.gridDegrees)
        do_day_night       = not options.keep_day_night_together
        do_multi_overpass  = options.allow_multiple_overpasses_per_cell

        temp_str = "will allow" if do_multi_overpass else "will not allow"
        LOG.debug("Space griding " + temp_str + " multiple overpasses per grid cell.")

        # determine the grid size in number of elements
        grid_lon_size      = int(math.ceil(360.0 / grid_degrees))
        grid_lat_size      = int(math.ceil(180.0 / grid_degrees))
        space_grid_shape   = (grid_lat_size, grid_lon_size) # I've confirmed with Nadia that this is the correct order
        
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
        expected_space_file_suffixes = io_manager.get_list_of_suffixes(DAILY_SPACE_TYPE, ALL_FILES_TYPE)
        for var_name in all_vars :

            for suffix in expected_space_file_suffixes :
                temp_stem = io_manager.build_name_stem(var_name, date_time=date_time_temp, satellite=satellite, suffix=suffix)
                temp_name = fbf.filename(temp_stem, TEMP_DATA_TYPE, shape=space_grid_shape)
                if os.path.exists(os.path.join(output_path, temp_name)) :
                    LOG.warn ("Cannot process files because matching temporary or output files exist in the output directory.")
                    return

        # loop to deal with data from each of the files
        failed_files       = 0
        successful_files   = 0
        abstract_data_sets = io_manager.get_expected_abstract_sets(instrument, separate_day_night=do_day_night)
        collected_data     = { }
        for each_file in sorted(possible_files) :

            full_file_path = os.path.join(input_path, each_file)

            LOG.debug("Processing file: " + full_file_path)

            # load the aux data
            file_object, temp_aux_data = io_manager.load_aux_data(full_file_path, min_scan_angle)
            # figure out what data sets we need to process
            data_sets = io_manager.get_expected_data_sets_from_aux_data (instrument, temp_aux_data, do_separate_day_night=do_day_night)

            ok_file     = True
            lon_indices = { }
            lat_indices = { }
            try :

                # calculate the indices for the space grid based on the navigation data
                # (we can do this now since the lon/lat is the same for each variable in the file)
                for set_key in data_sets.keys() :

                    set_mask      = data_sets[set_key][SET_MASK_KEY]
                    temp_lon_data = data_sets[set_key][LON_KEY][set_mask]
                    temp_lat_data = data_sets[set_key][LAT_KEY][set_mask]

                    lat_index, lon_index = space_gridding.calculate_index_from_nav_data(temp_lat_data, temp_lon_data, grid_degrees)
                    lat_indices[set_key] = lat_index
                    lon_indices[set_key] = lon_index

            except Exception, e :

                LOG.warn("Unable to process basic space griding for file: " + full_file_path)
                LOG.warn("This file will not be processed.")

                exc_type, exc_value, exc_traceback = sys.exc_info()
                LOG.debug(traceback.format_exception(exc_type, exc_value, exc_traceback))

                ok_file       = False
                failed_files += 1

            # if the file looks alright so far, continue processing it
            if ok_file :

                # loop to load each variable in the file and process it
                for variable_name in expected_vars[each_file] :

                    LOG.debug("Processing variable: " + variable_name)

                    # load the variable
                    file_object, var_data = io_manager.load_variable_from_file (variable_name,
                                                                                file_path=full_file_path,
                                                                                file_object=file_object)

                    # split the variable data by sets
                    separated_data   = { }
                    separated_time   = { }
                    separated_angles = { }
                    for set_key in data_sets.keys() :

                        separated_data  [set_key] =                                    var_data[data_sets[set_key][SET_MASK_KEY]]
                        separated_time  [set_key] = data_sets[set_key][SCAN_LINE_TIME_KEY]     [data_sets[set_key][SET_MASK_KEY]]
                        separated_angles[set_key] = data_sets[set_key][SENSOR_ZENITH_ANGLE_KEY][data_sets[set_key][SET_MASK_KEY]]

                    ok_file = True
                    space_grids  = { }
                    density_maps = { }
                    nobs         = { }
                    max_depths   = { }
                    aux_times    = { }
                    aux_angles   = { }
                    try :

                        # space grid the data using the indexes we calculated earlier
                        for set_key in data_sets.keys() :

                            # note: also preserve useful aux data for this file
                            temp_space_grid, temp_density_map, temp_nobs, temp_max_depth, temp_aux_time, temp_aux_angle = \
                                                                space_gridding.space_grid_data(grid_lat_size,
                                                                                               grid_lon_size,
                                                                                               separated_data[set_key],
                                                                                               lat_indices[set_key],
                                                                                               lon_indices[set_key],
                                                                                               aux_time=separated_time[set_key],
                                                                                               aux_sensor_zenith_angle=separated_angles[set_key])
                            space_grids [set_key] = temp_space_grid
                            density_maps[set_key] = temp_density_map
                            nobs        [set_key] = temp_nobs
                            max_depths  [set_key] = temp_max_depth
                            aux_times   [set_key] = temp_aux_time  # save the avg time for each cell
                            aux_angles  [set_key] = temp_aux_angle # save the max angle for each cell

                    except Exception, e :

                        LOG.warn("Unable to process variable data space griding for file: " + full_file_path)
                        LOG.warn("This variable will not be processed.")

                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        LOG.debug(traceback.format_exception(exc_type, exc_value, exc_traceback))

                        ok_file       = False
                        failed_files += 1

                    # if the data in the file looks ok so far, save it to the output
                    if ok_file :

                        #print("space grid shape: " + str(space_grids[set_key].shape))

                        # save the space grids and density info for this variable and it's density map to files
                        for set_key in data_sets.keys() :

                            if do_multi_overpass :

                                # save temporary data to accumulate it as we go through all the files for a day

                                # save the gridded data
                                io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                        satellite=satellite,
                                                                                        suffix=set_key + "-" + TEMP_SUFFIX_KEY),
                                                            space_grid_shape, output_path, space_grids[set_key], TEMP_DATA_TYPE)
                                # save the grid density map
                                io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                        satellite=satellite,
                                                                                        suffix=set_key + "-" + DENSITY_SUFFIX + "-" + TEMP_SUFFIX_KEY),
                                                            space_grid_shape, output_path, density_maps[set_key], TEMP_DATA_TYPE)
                                # save the number of observations grid
                                io_manager.save_data_to_file(io_manager.build_name_stem (variable_name, date_time=date_time_temp,
                                                                                        satellite=satellite,
                                                                                        suffix=set_key + "-" + NOBS_SUFFIX + "-" + TEMP_SUFFIX_KEY),
                                                            space_grid_shape, output_path, nobs[set_key], TEMP_DATA_TYPE)

                            else :

                                # if we haven't processed this variable yet, add a dictionary for it
                                if variable_name not in collected_data :
                                    collected_data[variable_name] = { }

                                # if we have no measurements, expand the array to depth 1 to make numpy happy
                                space_grids[set_key] = _expand_array_if_needed(space_grids[set_key], 1)

                                # if there isn't any data for this set in our collection, just put what we have in to start with
                                if set_key not in collected_data[variable_name] :

                                    collected_data[variable_name][set_key] = { }
                                    current_set = collected_data[variable_name][set_key]

                                    # save the the 2D arrays
                                    current_set["density"] = density_maps[set_key]
                                    current_set["nobs"]    = nobs[set_key]
                                    current_set["times"]   = aux_times[set_key]
                                    current_set["angles"]  = aux_angles[set_key]

                                    # save the space gridded data (the 3D array)
                                    new_depth = int(space_grids[set_key].shape[0] * ARRAY_GROWTH_FACTOR) # expand the arrays a little extra
                                    current_set["space-gridded-data"] = _expand_array_if_needed(space_grids[set_key], new_depth)

                                else : # when we already have data for this set key, incorporate the new overpass appropriately

                                    # there are several possible cases:
                                    #           we have no data in that cell of the grid           <- use data from the new file
                                    #           we have data in that cell, and it's the same orbit <- add the new data to the end of the old data
                                    #           we have data in that cell, it's a diff orbit       <- either replace the data in that cell or ignore the new data
                                    #                    (whether you replace or ignore depends on whether the new or old data has the worst sensor zenith angle)

                                    # for convenience
                                    current_set = collected_data[variable_name][set_key]

                                    # pre-calculate where there is any data at all in our old and new data sets
                                    have_old_data    = current_set["nobs"] > 0
                                    have_new_data    = nobs[set_key]       > 0
                                    both_have_data   = have_old_data & have_new_data

                                    # some other calculations to support our masking
                                    better_new_angle = aux_angles[set_key] < current_set["angles"]
                                    time_diff        = numpy.abs(current_set["times"] - aux_times[set_key])

                                    # figure out the masks that will control how we change our data

                                    # mask of the places where there is data in the new file, but we had none before
                                    have_only_new_data_mask = (~ have_old_data) & have_new_data
                                    # mask of the places where there is data in both and it's the same orbits
                                    use_both_mask = both_have_data & (time_diff <= space_gridding.SAME_TIME_RANGE_SECONDS)
                                    # mask of the places where there is data in both and it's different orbit
                                    use_only_new_data_mask = both_have_data & (time_diff > space_gridding.SAME_TIME_RANGE_SECONDS) & better_new_angle
                                    # Note: We will choose the orbit with the smallest maximum observed sensor zenith angle
                                    # in the grid cell – especially necessary at high latitudes
                                    use_new = have_only_new_data_mask | use_only_new_data_mask

                                    # expand the arrays if needed
                                    o_depth   = current_set["space-gridded-data"].shape[0] # the depth of the old array
                                    n_depth   = space_grids[set_key].shape[0]              # the depth of the new array
                                    c_depth   = numpy.max(current_set["density"][use_both_mask] + density_maps[set_key][use_both_mask]) if numpy.any(use_both_mask) else 0 # the combined depth
                                    new_depth = o_depth   if o_depth   >= n_depth else int(n_depth * ARRAY_GROWTH_FACTOR)
                                    new_depth = new_depth if new_depth >= c_depth else int(c_depth * ARRAY_GROWTH_FACTOR)
                                    new_space                         = _expand_array_if_needed(space_grids[set_key],              new_depth)
                                    current_set["space-gridded-data"] = _expand_array_if_needed(current_set["space-gridded-data"], new_depth)

                                    # replace any data where we are going to use just the new set
                                    current_set["space-gridded-data"][:, use_new] =          new_space[:, use_new]
                                    current_set["times"]                [use_new] =    aux_times[set_key][use_new]
                                    current_set["angles"]               [use_new] =   aux_angles[set_key][use_new]
                                    current_set["density"]              [use_new] = density_maps[set_key][use_new]
                                    current_set["nobs"]                 [use_new] =         nobs[set_key][use_new]

                                    # combine the data where we want to use both sets TODO, how can I do this in a more numpy and python friendly way?
                                    temp_shape = current_set["space-gridded-data"].shape
                                    for lat in range(temp_shape[1]) :
                                        for lon in range(temp_shape[2]) :
                                            if use_both_mask[lat, lon] :
                                                prev_num     = current_set["density"][lat, lon]
                                                num_adding   = density_maps[set_key] [lat, lon]
                                                new_total    = prev_num + num_adding
                                                current_set["space-gridded-data"][prev_num:new_total, lat, lon] =  new_space[:num_adding, lat, lon]
                                    # build a combined average of the times
                                    current_set["times"]  [use_both_mask] =  ((aux_times[set_key][use_both_mask] * nobs[set_key][use_both_mask]) + \
                                                                             (current_set["times"][use_both_mask] * current_set["nobs"][use_both_mask])) \
                                                                             / (nobs[set_key][use_both_mask] + current_set["nobs"][use_both_mask])
                                    temp_new_angles = aux_angles[set_key][use_both_mask]
                                    temp_old_angles = current_set["angles"][use_both_mask]
                                    current_set["angles"] [use_both_mask] = numpy.where(temp_new_angles > temp_old_angles, temp_new_angles, temp_old_angles) # select the largest angle from the two sets
                                    current_set["density"][use_both_mask] += density_maps[set_key][use_both_mask]
                                    current_set["nobs"]   [use_both_mask] +=         nobs[set_key][use_both_mask]

                        # if we got to here we processed the file successfully
                        successful_files += 1

            # make sure each file is closed when we're done with it
            io_manager.close_file(full_file_path, file_object)

        LOG.debug("Successfully processed " + str(successful_files) + " files and failed to process " + str(failed_files) + " files for this day.")
        
        # warn the user if we have fewer files than we need for this instrument
        if successful_files < (expected_num_files * EXPECTED_FRACTION_OF_FILES_PER_DAY) :
            LOG.warn("Processed " + str(successful_files)    + " files successfully for this day.")
            LOG.warn("Expected  " + str(expected_num_files) + " files for this instrument type.")
            
            if options.overrideMinCheck :
                LOG.warn ("Daily file(s) will be produced, but data may be unusable for this day.")
            else :
                LOG.critical("Daily file(s) will not be produced for this day due to lack of data." +
                             "If you wish to produce the daily file(s), rerun the program using the \'-p\' option.")
        
        # only save the daily data if we have enough files or have turned off the minimum check
        if ( (successful_files >= (expected_num_files * EXPECTED_FRACTION_OF_FILES_PER_DAY)) or options.overrideMinCheck ):

            # open the output workspace to load temp files if needed
            var_workspace = Workspace.Workspace(dir=output_path)

            # for use in labeling our output
            date_time_str = io_manager.get_datestamp (date_time_temp)

            # create a netCDF file to put our output in
            file_title = date_time_str + "_" + satellite + "_" + DAILY_SPACE_TYPE
            nc_output_file = io_manager.create_netCDF_output_file (output_path, file_title, do_ovewrite=False)

            # set up the file's default dimensions
            io_manager.set_up_dimensions_and_global_attrs_in_netCDF (nc_output_file,
                                                                     grid_lat_size, grid_lon_size,
                                                                     DAILY_SPACE_TYPE,
                                                                     date_time_str,
                                                                     global_attrs={"Satellite": satellite})

            # process each variable
            for variable_name in all_vars :

                if do_multi_overpass :

                    # collapse the per variable space grids to remove excess NaNs
                    LOG.debug("Packing space data for variable: " + variable_name)

                    for set_key in abstract_data_sets.keys( ) :

                        LOG.debug("Packing data for set key: " + set_key)

                        # load the density
                        temp_density = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                 satellite=satellite,
                                                                                 suffix=set_key + "-" + DENSITY_SUFFIX + "-" + TEMP_SUFFIX_KEY)][:]

                        # only process the final data if it exists
                        if numpy.sum(temp_density) > 0 :

                            # load the sparse space grid
                            var_data = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                satellite=satellite,
                                                                                suffix=set_key + "-" + TEMP_SUFFIX_KEY)][:]

                            # collapse the space grid
                            final_data = space_gridding.pack_space_grid(var_data, temp_density)

                            # save the final collapsed data to our netCDF output file
                            io_manager.add_variable_to_netCDF (nc_output_file,
                                                               variable_name,
                                                               final_data,
                                                               suffix_list = [set_key,])

                            # load the nobs file
                            nobs_counts = var_workspace[io_manager.build_name_stem(variable_name, date_time=date_time_temp,
                                                                                         satellite=satellite,
                                                                                         suffix=set_key + "-" + NOBS_SUFFIX + "-" + TEMP_SUFFIX_KEY)][:]

                            # collapse the nobs for the whole day
                            nobs_final = numpy.sum(nobs_counts, axis=0)

                            # save the final nobs to our netCDF output file
                            io_manager.add_variable_to_netCDF (nc_output_file,
                                                               variable_name,
                                                               nobs_final,
                                                               suffix_list = [set_key, NOBS_SUFFIX])

                        else :
                            LOG.warn("No " + set_key + " data was found for variable " + variable_name + ". Corresponding files will not be written.")

                else : # when we are limiting to a single overpass per cell

                    # save the collected space gridded data to files
                    LOG.debug("Saving space gridded data for variable: " + variable_name)

                    for set_key in collected_data[variable_name] :

                        LOG.debug("Saving data for set key: " + set_key)

                        # for convenience
                        current_set = collected_data[variable_name][set_key]

                        # figure out the max depth for the 3D files
                        max_depth = numpy.nanmax(current_set["density"])

                        # only process the final data if it exists
                        if max_depth > 0 :

                            # save the space gridded data to our netCDF output file
                            io_manager.add_variable_to_netCDF (nc_output_file,
                                                               variable_name,
                                                               current_set["space-gridded-data"][0:max_depth, :, :],
                                                               suffix_list = [set_key,])

                            # save the final nobs to our netCDF output file
                            io_manager.add_variable_to_netCDF (nc_output_file,
                                                               variable_name,
                                                               current_set["nobs"],
                                                               suffix_list = [set_key, NOBS_SUFFIX])

                        else :
                            LOG.warn("No " + set_key + " data was found for variable " + variable_name + ". Corresponding files will not be written.")

        # remove the extra temporary files in the output directory
        temp_suffix_list = io_manager.get_list_of_suffixes(DAILY_SPACE_TYPE, TEMP_FILE_TYPE)
        remove_suffixes = ["*" + p + "*" for p in temp_suffix_list]
        _remove_file_patterns(output_path, remove_suffixes)
    
    def time_griding_day(*args) :
        """given a day worth of files of daily space gridded data, calculate daily stats
        
        given an input directory that contains space gridded files for a day,
        calculate daily stats and put the resulting time gridded files
        for that day in the output directory.
        """
        
        # set up some of our caller determined settings for easy access
        input_path      = stg_util.clean_path(options.inputPath)
        output_path     = stg_util.clean_path(options.outputPath)
        stg_util.setup_dir_if_needed(output_path, "output")
        fix_nobs_cutoff = options.fixedNobsCutoff   if options.fixedNobsCutoff   >= 0 else None
        dyn_nobs_cutoff = options.dynamicNobsCutoff if options.dynamicNobsCutoff >= 0 else None
        nobs_LUT_path   = options.nobsLUT

        # get the nobs LUT if it was provided
        # TODO, this isn't properly tested because I haven't made a test nobs look up table yet
        # TODO, this also needs to move to a netCDF input
        nobs_LUT = None
        if nobs_LUT_path is not None :
            path_temp = os.path.split(nobs_LUT_path)
            temp_workspace = Workspace.Workspace(dir=path_temp[0])
            file_stem_temp = path_temp[0].split(".")[0]

            nobs_LUT = temp_workspace[file_stem_temp][:]

        # get a rough list of the daily gridded files
        possible_files = [ ]
        for file_name in os.listdir(input_path) :
            print("file name: " + file_name)
            if io_manager.is_stg_output_file(file_name, specific_type=DAILY_SPACE_TYPE) :
                possible_files.append(os.path.join(input_path, file_name))

        # process each daily gridded file
        for possible_file in possible_files :

            LOG.debug("Processing daily space gridded file: " + possible_file)

            # open the file object
            file_object = Dataset(possible_file, mode="r") # TODO, do I need to catch when it fails to open files?

            # pull some general information from the file
            date_time_str = file_object.DateTime
            satellite     = file_object.Satellite
            lat_size, lon_size = io_manager.get_nc_latlon_sizes(file_object)

            # sort the variables so they are organized by name and categories
            sorted_var_names = io_manager.sort_variable_names(file_object)
            # sorted format: {var name: {"cats": <category list string>, "nobs": the nobs variable name}}
            #print("variables: " + str(file_object.variables.keys()))
            #print("sorted:    " + str(sorted_var_names))

            # create our output file
            file_title = date_time_str + "_" + satellite + "_" + DAILY_TIME_TYPE
            nc_output_file = io_manager.create_netCDF_output_file (output_path, file_title, do_ovewrite=False)

            # set up the file's default dimensions
            io_manager.set_up_dimensions_and_global_attrs_in_netCDF (nc_output_file,
                                                                     lat_size, lon_size,
                                                                     DAILY_TIME_TYPE,
                                                                     date_time_str,
                                                                     global_attrs={"Satellite": satellite})

            # process the data for each variable
            for var_name in sorted_var_names.keys() :

                LOG.debug("Processing variable: " + var_name)

                # get the variable object for convenience
                var_object = file_object.variables[var_name]

                # load some info on the variable
                orig_name = var_object.originalVarName

                # load the raw data and the nobs
                gridded_data = var_object[:]
                nobs_data    = file_object.variables[sorted_var_names[var_name]["nobs"]][:]

                # build the cutoff mask
                bad_data     = time_gridding.create_sample_size_cutoff_mask(nobs_data,
                                                                            fixed_cutoff=fix_nobs_cutoff,
                                                                            dynamic_std_cutoff=dyn_nobs_cutoff,
                                                                            nobs_lut=nobs_LUT,)
                # apply the cutoff mask
                clean_gridded_data           = gridded_data.copy()
                clean_gridded_data[bad_data] = numpy.nan

                # figure out if we need to split our data
                masks_to_split = general_guidebook.mask_variable_for_time_gridding(satellite, orig_name, clean_gridded_data)

                # for each mask given, analyze the data selected by that mask
                for mask_key in masks_to_split.keys() :

                    if len(mask_key) > 0 :
                        LOG.debug("Processing masked subset: " + mask_key)

                    # select only the data from this mask, with the rest set to be nan
                    this_mask                 = masks_to_split[mask_key]
                    this_mask_data            = numpy.ones(clean_gridded_data.shape, dtype=TEMP_DATA_TYPE) * numpy.nan
                    this_mask_data[this_mask] = clean_gridded_data[this_mask]

                    # calculate the data fraction
                    num_mes          = numpy.sum(numpy.isfinite(this_mask_data), axis=0)
                    nobs             = nobs_data

                    # calculate the std
                    std_values  = numpy.nanstd(this_mask_data, axis=0)
                    # calculate the mean (divide by the number of measurements)
                    mean_values = numpy.nansum(this_mask_data, axis=0) / num_mes
                    # calculate the cloud fraction (num measurements / num obs)
                    cloud_frac  = num_mes / nobs
                    # calculate the uncertainty (std / num meas)
                    uncertainty = std_values / num_mes

                    # save the various stats to files

                    # make the new category list by adding the mask info to it
                    category_list = sorted_var_names[var_name]["cats"].split(" ")
                    if len(mask_key) > 0 :
                        category_list.append(mask_key)

                    # save the std and the mean
                    std_cats = list(category_list)
                    std_cats.append(STD_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      std_values,
                                                      suffix_list=std_cats)
                    mean_cats = list(category_list)
                    mean_cats.append(MEAN_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      mean_values,
                                                      suffix_list=mean_cats)

                    # save the cloud fraction and uncertainty
                    cloud_frac_cats = list(category_list)
                    cloud_frac_cats.append(CLOUD_FRACTION_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      cloud_frac,
                                                      suffix_list=cloud_frac_cats)
                    uncertainty_cats = list(category_list)
                    uncertainty_cats.append(UNCERTAINTY_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      uncertainty,
                                                      suffix_list=uncertainty_cats)

                    # save the number of measurements and observations
                    num_mes_cats = list(category_list)
                    num_mes_cats.append(NUM_MES_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      num_mes,
                                                      suffix_list=num_mes_cats)
                    nobs_cats = list(category_list)
                    nobs_cats.append(NOBS_SUFFIX)
                    io_manager.add_variable_to_netCDF(nc_output_file,
                                                      orig_name,
                                                      nobs,
                                                      suffix_list=nobs_cats)

    def time_griding_multiday(*args) :
        """given a directory with multiple days of daily space gridded data, calculate overall stats

        given an input directory that contains appropriate daily stats for
        more than one day, calculate overall stats for that time period and
        put the resulting gridded files for that time period in the output
        directory.

        Note: The input days are expected to be consecutive (possibly with missing days).
        The program will not check the dates present. Sparse (in time) data may give
        unexpected results.
        """

        # TODO, this is unfinished

    def make_nobs_lut (*args) :
        """given a directory with a multiple daily space gridded files, make a nobs look up table

        generally this will expect a month of daily space gridded files and will
        output some files with statistical information on the nobs across that period
        the resulting look up tables are intended to be used with dynamic cutoffs in
        time_griding_multiday calls
        """

        # TODO, transition this to reading in netCDF daily files instead

        # set up some of our caller determined settings for easy access
        input_path      = stg_util.clean_path(options.inputPath)
        output_path     = stg_util.clean_path(options.outputPath)
        stg_util.setup_dir_if_needed(output_path, "output")

        # check the directory for sets of daily files
        expected_files_by_date = defaultdict(list)
        possible_files = os.listdir(input_path)
        for file_name in sorted(possible_files) :
            if file_name.find(NOBS_SUFFIX) < 0 :
                LOG.debug("Disregarding non-nobs file: " + str(file_name))
            else :
                date_stamp = io_manager.get_date_stamp_from_file_name(file_name)
                if date_stamp is None :
                    LOG.debug("Disregarding file with no date stamp: " + str(file_name))
                else :
                    expected_files_by_date[date_stamp].append(file_name)

        # organize the daily files
        organized_files = { }
        for date_stamp in expected_files_by_date.keys() :

            organized_files[date_stamp] = io_manager.organize_space_gridded_files(expected_files_by_date[date_stamp])

        # set up the variable workspace so we can load our input files
        var_workspace = Workspace.Workspace(dir=input_path)

        # for each set of daily files
        for date_stamp in organized_files.keys() :
            for set_key in organized_files[date_stamp].keys() :

                # load the nobs
                nobs_stem    = organized_files[date_stamp][set_key][NOBS_SUFFIX].split('.')[0]
                nobs_data    = var_workspace[nobs_stem][:]
                stem_no_time =    nobs_stem[nobs_stem.find("_")+1:] # get rid of the time
                stem_no_time = stem_no_time[0:stem_no_time.rfind("_")] # get rid of the previous suffix

                # save the number of observations grid
                io_manager.save_data_to_file(stem_no_time + "_" + set_key + "-" + NOBS_LUT_SUFFIX,
                                            nobs_data[0].shape, output_path, nobs_data[0], TEMP_DATA_TYPE)
                io_manager.save_data_to_file(stem_no_time + "_" + ALL_SET_KEY + "-" + NOBS_LUT_SUFFIX,
                                            nobs_data[0].shape, output_path, nobs_data[0], TEMP_DATA_TYPE)
                # TODO, does the all set key need to be collapsed? also, what if stuff is already in the all set?

    ##### This is the end of the menu selection functions. #####

    # all the local public functions are considered part of the application, collect them up
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    # if what the user asked for is not one of our existing functions, print the help
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        # help()
        return 9
    else:
        start_time = time.time()

        # call the function the user named, given the arguments from the command line
        rc = locals()[args[0]](*args[1:])

        end_time = time.time()

        LOG.info("Elapsed seconds during this run: " + str(end_time - start_time))

        return 0 if rc is None else rc

if __name__=='__main__':
    sys.exit(main())

