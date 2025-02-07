
import time
import glob
import gzip
import re
import os

import uuid

import parse_embl

###############################################################################
# Functions used as Dask Tasks
###############################################################################

def glob_subdirs(dir_path: str) -> tuple:
    """
    search for subdirectories in the provided directory path

    Parameters
    ----------
        dir_path
            string, global or local path within which the search for subdirs
            will occur. 

    Returns
    -------
        "glob_subdirs"
            string used to ID type of task
        subdirs
            list of strings corresponding to the found subdirs
        `time.time() - st`
            elapsed time for this task, units: seconds
        dir_path
            same as given input
    """
    st = time.time()
    return "glob_subdirs", glob.glob(dir_path + "/*/"), time.time() - st, dir_path


def glob_files(dir_path: str) -> tuple:
    """
    return list of files matching the search string
    
    Parameters
    ----------
        dir_path
            string, global or local path within which the search for subdirs
            will occur.

    Returns
    -------
        "glob_files"
            string used to ID type of task
        files
            list of strings corresponding to the found files
        `time.time() - st`
            elapsed time for this task, units: seconds
        dir_path
            same as given input
    """
    st = time.time()
    # grab all file path strings in the given dir_path
    files = glob.glob(dir_path + f"/*dat.gz")
    
    # only a subset of data files in the ENA sequence/ subdir are of interest 
    # to us.
    # THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    if "sequence" in dir_path:
        pattern = re.compile(r".*(ENV|PRO|FUN|PHG).*")
        # NOTE: regex to only gather file names with (ENV|PRO|FUN|PHG)
        files = [file_ for file_ in files if pattern.match(file_)]

    return "glob_files", files, time.time() - st, dir_path
    #return [(file_name, os.path.getsize(file_name)) for file_name in glob.glob(dir_path + f"/{search_string}")]


def process_many_files(
        file_path_list: list, 
        database_params: dict,
        db_name: str, 
        common_output_dir: str):
    """
    given a list of files, process them one at a time. Gather the files written 
    during processing and return that list. 

    PARAMETERS
    ----------
        file_path_list
            list of strings or pathlib.Path objs, assumed to be associated with
            gzipped EMBL/GenBank flat files. 
        database_params 
            dict or configparser.ConfigParser obj. necessary keys or 
            attributes are "user", "password", "host", "port"
        db_name
            string, name of the EFI database to be used to perform
            queries. 
        common_output_dir
            string, local or global path within which result files will be 
            written.

    RETURNS
    -------
        "process_many_files"
            string used to ID type of task
        tab_files
            list of strings corresponding to the tsv files written during the
            tawk
        `time.time() - st`
            elapsed time for this task, units: seconds
        file_path_list
            same as given input

    """
    st = time.time()
    task_uuid = uuid.uuid4()
    #task_uuid = uuid.uuid1()
    out_dir = common_output_dir + f"/{task_uuid}"
    os.makedirs(out_dir)
   
    # connect to the database
    db_connection = IDMapper(database_params,db_name)

    tab_files = []
    for file_ in file_path_list:
        start_time = time.time()
        # instead use regex to parse the file name so that the path of parsed
        # file can be recreated for the output file name
        fn_name = Path(file_).name.split('.')[0]
        tab_file = parse_embl.process_file(file_, db_connection, out_dir + f"/{fn_name}.tab")
        stop_time = time.time()
        # if the file_ does not return tab results, no file will be written
        if os.path.isfile(out_dir + f"/{fn_name}.tab"):
            tab_files.append(tab_file)
    
    db_connection.close()
    
    # quick and dirty handling of the task output
    if not tab_files:
        tab_files = [""]
    
    return "process_many_files", tab_files, time.time() - st, file_path_list


