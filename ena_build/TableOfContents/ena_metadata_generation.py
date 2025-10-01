
import os
import shutil
import sys
import argparse
import time
import json

import dask
from distributed import Client, as_completed

import mysql_database
from ..workflow_logging import setup_logger, clean_logger
from ..glob_tasks import glob_subdirs, glob_files
from toc_tasks import gather_files_metadata

###############################################################################
# Parse Input Arguments and Files
###############################################################################

def parse_input_arguments() -> argparse.Namespace:
    """
    Returns an `<argparse.Namespace>` with attributes associated with the input
    arguments for the ENA database build script:
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
        * ``--md5-hash`` or ``-md5``, flag to calculate the md5 hash for each
                                     file.
        * ``--id-index`` or ``-index``, flag to control whether assembly files
                                        are processed to gather protein_ids
                                        contained in the associated file.

    Returns
    -------
        An :external+python:py:class:`argparse.Namespace` object with attributes
            args.ena_paths
            args.output_dir
            args.scheduler_file
            args.n_workers
            args.tskmgr_log_file
            args.md5_hash
            args.id_index
    """
    parser = argparse.ArgumentParser(
        description = "Create a Table of Contents for the ENA Database"
    )
    parser.add_argument("--ena-paths", required=True, nargs="+", help="Arbitrary number of file paths that house subdirectories to be searched for ENA related dat.gz files.")
    parser.add_argument("--output-dir", "-out", required=True, help="Path to the common output directory within which subdirectories and associated tab-separated data files will be saved.")
    parser.add_argument("--scheduler-file", "-s", help="Path string to the dask scheduler file.")
    parser.add_argument("--n-workers", "-nWorkers", default = 2, type=int, help="Number of workers available to perform tasks, default = 2.")
    parser.add_argument("--tskmgr-log-file", "-log", default = "dask_tskmgr.log", help="Path string for a logging file, default = 'dask_tskmgr.log'.")
    parser.add_argument("--md5-hash", "-md5", action = "store_true", help="Flag to set whether the md5 hash is calculated for each file in the TOC.")
    parser.add_argument("--id-index", "-index", action = "store_true", help="Flag to control whether assembly files are processed to gather protein_ids contained in files.")
    args = parser.parse_args()
    return args


###############################################################################
# WORKFLOW FUNCTION
###############################################################################

def workflow():
    """ Run the Dask workflow to parse the ENA dataset. """

    # parse input arguments
    args = parse_input_arguments()

    # make the args.output_dir
    os.makedirs(args.output_dir, exist_ok=True)

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

    processing_tasks = 0
    finished_processing_tasks = 0

    # NOTE: this is a rough optimization for handling large lists of unevenly
    # distributed (in subdirectories) sets of files.
    ideal_nFiles = 10

    ## pre-compile the regex patterns
    # source_pattern is used as a file filter to only consider files from the
    # given sources
    source_pattern = re.compile(r"_(ENV|PRO|FUN|PHG)_")
    # dir_pattern is used to parse subdirectories' names to make reporting in
    # TOC/id index agnostic to absolute file paths
    dir_pattern = re.compile(r"(wgs)\/(\w*)\/(\w*)|(sequence)\/(\w*)")
    # file_name_pattern is used to get the root name of the file
    file_name_pattern = re.compile(r"\/(\w*)\.dat\.gz")

    # submit tasks to the client that glob search for the intermediate layer
    # of subdirs in ENA directory tree
    glob_subdirs_futures = client.map(glob_subdirs, args.ena_paths)
    # setup the iterator that is filled with futures as they complete; this
    # tasks_completed object also lets us add new tasks to the queue, making
    # this for loop very flexible/dynamic.
    tasks_completed = as_completed(glob_subdirs_futures)
    # loop over finished tasks
    with open(args.output_dir + "ENA_file_toc.tab","w") as tab:
        for finished_task in tasks_completed:
            # gather the results from the finished task
            # results is the return tuple from any of the task functions, so
            # handle them accordingly.
            results = finished_task.result()
            if not results[1]:
                # results[1] will always be some true-equivalent value unless
                # if a glob search returns an empty list. if this is the case,
                # then move on. No new tasks need to be submitted. Log the
                # result.
                if "glob" in results[0]:
                    main_logger.info(
                        f"{results[-1]} did not have any expected files/subdirs"
                    )
                if results[0] == "gather_files_metadata":
                    finished_processing_tasks += 1
                    main_logger.info(
                        f"No file metadata was gathered for the set of files."
                        + f" This took {results[2]} seconds."
                        + f" {finished_processing_tasks} tasks completed out of"
                        + f" {processing_tasks}."
                    )
                else:
                    main_logger.info(
                        "Something's gone wrong. Closing down. Time: "
                        + f"{time.time()}"
                    )
                    clean_logger(main_logger)

            elif results[0] == "glob_subdirs":
                # finished task is a glob_subdirs task, so results[1] will be
                # a list of subdirectories. For each subdir, submit a new task
                # to glob for gzipped files.
                main_logger.info(
                    f"Found {len(results[1])} subdirectories in {results[3]}."
                    + f" Took {results[2]} seconds. Submitting "
                    + f"{len(results[1])} new tasks to search for gzipped "
                    + "files.")
                # list comprehension is submitting a new task to the client,
                # one for each subdirectory found in the intermediate
                # directory. The new future gets added to the task_completed
                # iterator so will be gathered and logged in this for loop.
                new_futures = [
                    client.submit(
                        glob_files,
                        subdir,
                        source_pattern
                    ) for subdir in results[1]
                ]
                for new_future in new_futures:
                    tasks_completed.add(new_future)

            elif results[0] == "glob_files":
                # finished task is a glob_files task, so results[1] will be
                # the list of gzipped files. Break this list down into bite
                # sized chunks and submit a new task for each chunk.
                #shards = [results[1][i::args.n_workers] for i in range(args.n_workers)]
                nTasks = int(len(results[1])/ideal_nFiles) + 1
                shards = [results[1][i::nTasks] for i in range(nTasks)]

                # list comprehension is submitting a new task to the client,
                # one for each worker, evenly separating the number of files
                # to be processed across the tasks. If n_workers is > than
                # files in results[1], then only the necessary number of tasks
                # to process one file per worker are created.
                #new_futures = [client.submit(process_many_files, shard, database_params = database_params, db_name = args.db_name, final_output_dir = args.output_dir, temp_output_dir = args.local_scratch) for shard in shards if shard]
                new_futures = [
                    client.submit(
                        gather_files_metadata,
                        shard,
                        toc_bool = True,
                        md5_hash_bool = args.md5_hash
                        index_bool = args.id_index
                    ) for shard in shards if shard
                ]
                main_logger.info(f"Found {len(results[1])} gzipped files in "
                    + f"{results[3]}. Took {results[2]} seconds. Sharding the "
                    + f"list into {len(new_futures)} tasks.")
                processing_tasks += len(new_futures)
                for new_future in new_futures:
                    tasks_completed.add(new_future)

            elif results[0] == "gather_files_metadata":
                # finished task is a gather_files_metadata task, so results[1]
                # is a list of FileMetadata objects, each containing metadata
                # about the files associated with the task's shard.
                finished_processing_tasks += 1
                main_logger.info(
                    f"Parsing {len(results[1])} file(s) took {results[3]} "
                    + f" seconds. {finished_processing_tasks} tasks completed"
                    + f" out of {processing_tasks}."
                )
                # loop over the FileaMetadata objects and write their
                # information out to file
                for metadata in results[1]:
                    # create a directory subtree that specifies where
                    # the assembly file is positioned in an assumed
                    # ENA directory tree
                    dir_subtree = os.path.join(
                        *[
                            elem for elem in dir_pattern.findall(
                                metadata.file_path
                            ) if elem
                        ]
                    )
                    # grab the file name from the file_path
                    file_name = file_name_pattern.findall(
                        metadata.file_path
                    )[0]

                    # make the values in the .toc dict into a tab separated str
                    toc_string = "\t".join(
                        [metadata.toc[key] for key in toc_column_list]
                    )
                    
                    # write the TOC 
                    tab.write(f"{dir_subtree}\t{file_name}\t{toc_string}\n")

                    # write the index to a file as well
                    if args.id_index:
                        with open(
                            args.output_dir + "protein_id_index.tab","a"
                        ) as index:
                            # write the protein_id's index info to file
                            index.write(
                                "\n".join(
                                    [
                                        f"{id}\t{file_name}\t{dir_subtree}"
                                        for id in metadata.ids
                                    ]
                                )
                            )
                            index.write("\n")

    main_logger.info(
        f"Closing dask pipeline and logging. Time: {time.time()}"
    )
    clean_logger(main_logger)


###############################################################################
# MAIN
###############################################################################

if __name__ == "__main__":
    workflow()

