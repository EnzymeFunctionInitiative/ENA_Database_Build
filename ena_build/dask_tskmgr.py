
import sys
import argparse
import logging
import configparser
import time

import dask
from distributed import Client, as_completed

#import build_ena_db
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
        * ``--db-type`` or ``-type``, string denoting type of database file to 
                                      be queried
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
    # EFI database input arguments
    parser.add_argument("--db-type", "-type", required=True, help="string denoting type of database file to be queried, expected values are 'mysql' or 'sqlite'")
    parser.add_argument("--db-config", "-conf", required=True, help="file path to the config file for database connection; assumed format is Windows INI")
    parser.add_argument("--db-name", "-dbn", required=True, help="name of the EFI database to query for retrieving IDs")
    # ENA file IO input arguments
    parser.add_argument("--ena-paths", required=True, nargs="+", help="arbitrary number of file paths that house subdirectories to be searched for ENA related dat.gz files.")
    parser.add_argument("--output-dir", "-out", required=True, help="path to the common output directory within which subdirectories and parsed data will be saved.")
    # dask specific input arguments
    parser.add_argument("--scheduler-file", "-s", required=True, help="path string to the dask scheduler file")
    parser.add_argument("--n-workers", "-nWorkers", required=True, type=int, help="number of workers available to perform tasks")
    parser.add_argument("--tskmgr-log-file", "-log", default = "dask_taskmngr.log", help="path string for a logging file")
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
# MAIN
###############################################################################

if __name__ == "__main__":
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
    
    # NOTE: validate other input arguments... 
    # ensure that any paths are global

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
    # need to include at the moment because this file is not in the pythonpath
    #client.upload_file("/home/n-z/rbdavid/Projects/ENA_building/build_ena_db.py")
    
    # prep a list to gather final tab files
    tab_files = []
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
            main_logger.info(f"{results[-1]} did not have any expected files/subdirs")
        elif results[0] == "glob_subdirs":
            # finished task is a glob_subdirs task, so results[1] will be the
            # list of subdirectories. For each subdir, submit a new task to 
            # glob for gzipped files. 
            main_logger.info(f"Found {len(results[1])} subdirectories in {results[3]}. Took {results[2]} seconds. Submitting {len(results[1])} new tasks to search for gzipped files.")
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
            shards = [results[1][i::args.n_workers] for i in range(args.n_workers)]
            # list comprehension is submitting a new task to the client, one 
            # for each worker, evenly separating the number of files to be 
            # processed across the tasks. If n_workers is > than files in 
            # results[1], then only the necessary number of tasks to process 
            # one file per worker are created.
            new_futures = [client.submit(process_many_files, shard, database_params = database_params, db_name = args.db_name, common_output_dir = args.output_dir) for shard in shards if shard]
            main_logger.info(f"Found {len(results[1])} gzipped files in {results[3]}. Took {results[2]} seconds. Sharding the list into {len(new_futures)} tasks.")
            for new_future in new_futures:
                tasks_completed.add(new_future)
        elif results[0] == "process_many_files":
            # NOTE: UPDATE THIS
            # finished task is a process_many_files task, so results[1] is a 
            # list of file paths to temp tab files that need to be concatenated.
            main_logger.info(f"Parsing {len(results[3])} files took {results[2]} seconds. Individual tab files written to {args.output_dir}.")
            tab_files += results[1]

        # NOTE: delete me 
        else:
            main_logger.info(results)
            sys.exit("Unexpected results from a task. Exiting.")

    # gather all results into a final file?
    # or 
    # prep a bash script that concatenates the tab files and creates a table or
    # add a column to the mysql database?
    #for tab_file in tab_files:
    #    shutil.move(tab_file, args.output_dir + f"/{Path(tab_file).name}")

    ## for now, just log where the tab files were written
    main_logger.info(f"Find tab files in subdirs w/in {args.output_dir}")
    
    main_logger.info(f"Closing dask pipeline and logging. Time: {time.time()}")
    clean_logger(main_logger)

