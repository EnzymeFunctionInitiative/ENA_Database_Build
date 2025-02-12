
import time
import glob
import gzip
import re
import os

import mysql_database
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
    # use regex to match the parent directories' names; three layers worth if
    # in `wgs` tree of ENA or two layers worth if in `sequences` tree. This
    # regex will match a file path string, creating a list of a tuple with len
    # 7. First four elements are associated with the wgs tree, the remaining
    # with the sequences tree. 
    # assume that all files in the file_path_list are sourced from the same
    # directory
    # THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    dir_pattern = re.compile(r"(wgs)\/(\w*)\/(\w*)|(sequences)\/(\w*)")
    # use regex to match the file name stem from the given file path; will 
    # create a list of len 1. 
    file_pattern= re.compile(r"\/(\w*)\.dat\.gz")

    # only grab groups that were successfully matched. 
    matches = [elem for elem in dir_pattern.findall(file_path_list[0])[0] if elem]
    # create an output_dir string that easily maps to the files being parsed
    if common_output_dir[-1] != "/":
        common_output_dir += "/"
    out_dir = common_output_dir + "-".join(matches)
    # make the directory
    os.makedirs(out_dir,exist_ok=True)
   
    # connect to the database
    db_connection = mysql_database.IDMapper(database_params,db_name)

    tab_files = []
    for file_path in file_path_list:
        start_time = time.time()
        # grab the stem of the file name to use in writing results
        fn_name = file_pattern.findall(file_path)[0]
        # process the file
        tab_file = parse_embl.process_file(
            file_path, 
            db_connection, 
            out_dir + f"/{fn_name}.tab"
        )
        stop_time = time.time()
        # if the file does not return any tab results, no file will be written
        if os.path.isfile(out_dir + f"/{fn_name}.tab"):
            tab_files.append(tab_file)
    
    db_connection.close()
    
    # quick and dirty handling of the task output
    if not tab_files:
        tab_files = [""]
    
    return "process_many_files", tab_files, time.time() - st, file_path_list


