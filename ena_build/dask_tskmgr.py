
import os
import shutil
import sys
import argparse
import logging
import configparser
import time

import dask
from distributed import Client, as_completed

import mysql_database
from dask_tasks import glob_subdirs, glob_files, process_many_files

###############################################################################
# Logging Functions
###############################################################################

def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


###############################################################################
# Parse Input Arguments and Files
###############################################################################

def parse_input_arguments() -> argparse.Namespace:
    """
    Returns an `<argparse.Namespace>` with attributes associated with the input
    arguments for the ENA database build script:
        * ``--db-config`` or ``-conf``, path string to the database config file
        * ``--db-name`` or ``-dbn``, string for the database name
        * ``--ena-paths``, a number of path strings; accepts multiple values 
                           and returns a list of the values
        * ``--output-dir`` or ``-out``, path string within which files will be
                                        written
        * ``--scheduler-file`` or ``-s``, path string to the dask-distributed
                                            scheduler's json file for tracking 
                                            workers
        * ``--n-workers`` or ``-nWorkers``, number of dask-distributed workers
                                            available to perform tasks
        * ``--tskmgr-log-file`` or ``-log``, path string to be used for a log
                                             file that tracks progress
        * ``--local-scratch`` or ```-scratch``, path string where temp files 
                                                will be written, default = ''
                                                indicating to not write temp 
                                                files to storage
    
    Returns
    -------
        An :external+python:py:class:`argparse.Namespace` object with attributes
            args.db_config
            args.db_name
            args.ena_paths
            args.output_dir
            args.scheduler_file
            args.n_workers
            args.tskmgr_log_file
            args.local_scratch
    """
    parser = argparse.ArgumentParser(description = "Process the ENA Database")
    parser.add_argument("--db-config", "-conf", required=True, help="file path to the config file containing  database connection information; assumed format is Windows INI.")
    parser.add_argument("--db-name", "-dbn", required=True, help="name of the EFI database to query for retrieving IDs.")
    parser.add_argument("--ena-paths", required=True, nargs="+", help="arbitrary number of file paths that house subdirectories to be searched for ENA related dat.gz files.")
    parser.add_argument("--output-dir", "-out", required=True, help="path to the common output directory within which subdirectories and associated tab-separated data files will be saved.")
    parser.add_argument("--scheduler-file", "-s", help="path string to the dask scheduler file.")
    parser.add_argument("--n-workers", "-nWorkers", default = 2, type=int, help="number of workers available to perform tasks, default = 2.")
    parser.add_argument("--tskmgr-log-file", "-log", default = "dask_tskmgr.log", help="path string for a logging file, default = 'dask_tskmgr.log'.")
    parser.add_argument("--local-scratch", "-scratch", default = "", help="path string where temp files will be written, default = '', indicating do not write temp files to storage.")
    args = parser.parse_args()
    return args


def parse_config(file_path):
    """ Parse the database config file and return parameters """
    config = configparser.ConfigParser()
    try:
        config.read(file_path)
    except configparser.Error as err:
        if err.errno == configparser.ParsingError:
            sys.exit(f"Provided {file_path} file is not correctedly formatted."
                    + " This file should follow Windows INI formating.")
        else:
            sys.exit(f"Parsing --db-config file {file_path} failed:\n{err}")
    return config


###############################################################################
# WORKFLOW FUNCTION
###############################################################################

def workflow():
    """ Run the Dask workflow to parse the ENA dataset. """
    
    # parse input arguments
    args = parse_input_arguments()

    # parse --db-config and check for essential parameters
    database_params = parse_config(args.db_config)["database"]
    for param in ["user", "password", "host", "port"]:
        if param not in database_params.keys():
            sys.exit(f"'{param}' is missing from the --db-config file.")

    # test db connection before starting the whole pipeline
    idmapper = mysql_database.IDMapper(database_params, args.db_name)
    if not idmapper:
        sys.exit("Failed to connect to the MySQL database.")
    idmapper.close()
    
    # make the args.output_dir and args.local_scratch directories
    os.makedirs(args.output_dir, exist_ok=True)
    # args.local_scratch may be an empty string so check for that
    if args.local_scratch:
        os.makedirs(args.local_scratch, exist_ok=True)
        # check if a file can be written to args.local_scratch space and then
        # shutil'd to args.output_dir
        with open(args.local_scratch + '/test.txt','w') as out:
            out.write("hello world")
        out = shutil.move(args.local_scratch + '/test.txt', args.output_dir)
        os.remove(out)

    # set up the main logger file and list all relevant parameters
    main_logger = setup_logger("tskmgr_logger",args.tskmgr_log_file)
    main_logger.info(f"Starting dask pipeline and setting up logging. Time: {time.time()}")
    for arg in vars(args):
        main_logger.info(f"{arg}: {getattr(args,arg)}")

    # also list the dask parameters; this is only included for thoroughness 
    # sake; we haven't messed with any of these parameters
    dask_parameter_string = "#"*80 +"\nDask parameters:\n"
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}':'{value}'"
    dask_parameter_string += "\n" + "#"*80
    main_logger.info(f"\n{dask_parameter_string}")

    # start the dask client.
    client = Client(scheduler_file=args.scheduler_file)
    
    # prep a list to gather final tab files
    tab_files = []
    # prep count variables to track progress through the processing_many_files
    # tasks
    processing_tasks = 0
    finished_processing_tasks = 0
    ideal_nFiles = 10
    
    # submit tasks to the client that glob search for the intermediate layer 
    # of subdirs in ENA directory tree
    glob_subdirs_futures = client.map(glob_subdirs, args.ena_paths)
    # setup the iterator that is filled with futures as they complete; this 
    # tasks_completed object also lets us add new tasks to the queue, making 
    # this for loop very flexible/dynamic. 
    tasks_completed = as_completed(glob_subdirs_futures)
    # loop over finished tasks
    for finished_task in tasks_completed:
        # gather the results from the finished task
        # results is the return tuple from any of the task functions, so handle 
        # them accordingly.
        results = finished_task.result()
        if not results[1]:
            # results[1] will always be some true-equivalent value unless if a
            # glob search returns an empty list. if this is the case, then move
            # on. No new tasks need to be submitted. Log the result.
            main_logger.info(f"{results[-1]} did not have any expected " 
                + "files/subdirs")
        elif results[0] == "glob_subdirs":
            # finished task is a glob_subdirs task, so results[1] will be the
            # list of subdirectories. For each subdir, submit a new task to 
            # glob for gzipped files. 
            main_logger.info(f"Found {len(results[1])} subdirectories in "
                + f"{results[3]}. Took {results[2]} seconds. Submitting " 
                + f"{len(results[1])} new tasks to search for gzipped files.")
            # list comprehension is submitting a new task to the client, one 
            # for each subdirectory found in the intermediate directory. The
            # new future gets added to the task_completed iterator so will be
            # gathered and logged in this for loop.
            new_futures = [client.submit(glob_files, subdir) for subdir in results[1]]
            for new_future in new_futures:
                tasks_completed.add(new_future)
        elif results[0] == "glob_files":
            # finished task is a glob_files task, so results[1] will be the 
            # list of gzipped files. Break this list down into bite sized 
            # chunks and submit a new task for each chunk. 
            #shards = [results[1][i::args.n_workers] for i in range(args.n_workers)]
            nTasks = int(len(results[1])/ideal_nFiles) + 1
            shards = [results[1][i::nTasks] for i in range(nTasks)]

            # list comprehension is submitting a new task to the client, one 
            # for each worker, evenly separating the number of files to be 
            # processed across the tasks. If n_workers is > than files in 
            # results[1], then only the necessary number of tasks to process 
            # one file per worker are created.
            new_futures = [client.submit(process_many_files, shard, database_params = database_params, db_name = args.db_name, final_output_dir = args.output_dir, temp_output_dir = args.local_scratch) for shard in shards if shard]
            main_logger.info(f"Found {len(results[1])} gzipped files in " 
                + f"{results[3]}. Took {results[2]} seconds. Sharding the " 
                + f"list into {len(new_futures)} tasks.")
            processing_tasks += len(new_futures)
            for new_future in new_futures:
                tasks_completed.add(new_future)
        elif results[0] == "process_many_files":
            # finished task is a process_many_files task, so results[1] is a 
            # list of file paths of tab files that need to be concatenated.
            finished_processing_tasks += 1
            task_tab_files = [tab for tab in results[1] if tab]
            main_logger.info(f"Parsing {len(results[3])} file(s) took " 
                + f"{results[2]} seconds with {len(task_tab_files)} tab "
                + f"file(s) written. {finished_processing_tasks} tasks " 
                + f"completed out of {processing_tasks}.")
            tab_files += task_tab_files

    main_logger.info("Find individual tab files in subdirs w/in" 
        + f" {args.output_dir} .")
    
    # sort the tab_files list... likely to be slow since so many tab files 
    # will be written
    tab_files.sort()
    # open a file to be filled with all contents of all of the tab files
    with open(args.output_dir + '/ena.tab','w') as out:
        # loop over the sorted tab files
        for tab in tab_files:
            # open the tab file and write it out to the out file
            with open(tab,'r') as tsv:
                shutil.copyfileobj(tsv, out)

    main_logger.info("Catenated all individual tabs into" 
        + f" {args.output_dir + '/ena.tab'}.")

    # need to clean up scratch space before closing down the workflow
    if args.local_scratch:
        main_logger.info("Cleaning the scratch space.") 
        dirs = os.listdir(args.local_scratch)
        for direc in dirs:
            try:
                os.rmdir(direc)
            except:
                pass

    main_logger.info(f"Closing dask pipeline and logging. Time: {time.time()}")
    clean_logger(main_logger)


###############################################################################
# MAIN
###############################################################################

if __name__ == "__main__":
    workflow()

