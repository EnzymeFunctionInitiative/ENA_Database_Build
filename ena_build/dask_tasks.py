
import time
import glob
import gzip
import re
import os
import shutil

import mysql_database
import parse_embl

###############################################################################
# Functions used as Dask Tasks
###############################################################################

def glob_subdirs(dir_path: str) -> tuple:
    """
    Search for subdirectories in the provided directory path.

    Parameters
    ----------
        dir_path
            str, global or local path within which the search for subdirs
            will occur. 

    Returns
    -------
        "glob_subdirs"
            str, used to ID type of task.
        subdir_list
            list of strs, each element corresponding to a found subdir.
        `time.time() - st`
            float, elapsed time for this task, units: seconds.
        dir_path
            str, same as given input.
    """
    st = time.time()
    # Grab all subdirectory path strings in the given dir_path
    subdir_list = [
        dir_path + "/" + dir_.name 
        for dir_ in os.scandir(dir_path) 
        if not dir_.name.startswith('.') 
        and dir_.is_dir()
    ]
    return "glob_subdirs", subdir_list, time.time() - st, dir_path


def glob_files(dir_path: str) -> tuple:
    """
    Return list of files matching the search string.
    
    Parameters
    ----------
        dir_path
            str, global or local path within which the search for subdirs will
            occur.

    Returns
    -------
        "glob_files"
            str, used to ID type of task.
        files
            list of strs, each element corresponding to a found file.
        `time.time() - st`
            float, elapsed time for this task, units: seconds.
        dir_path
            str, same as given input.
    """
    st = time.time()
    # Grab all file path strings in the given dir_path
    files = [
        dir_path + "/" + file.name 
        for file in os.scandir(dir_path) 
        if file.name.endswith('.dat.gz') 
        and file.is_file()
    ]
    
    # Only a subset of data files in the ENA sequence/ subdir are of interest 
    # to us. As far as I know, the second underscored section of the file name
    # denote the origin species type, which is what we need to consider.
    # NOTE: THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    if "sequence" in dir_path:
        # NOTE: regex to only gather file names with (ENV|PRO|FUN|PHG) in them
        pattern = re.compile(r"_(ENV|PRO|FUN|PHG)_")
        files = [file_ for file_ in files if pattern.search(file_)]

    return "glob_files", files, time.time() - st, dir_path


def process_many_files(
        file_path_list: list, 
        database_params: dict,
        db_name: str, 
        final_output_dir: str,
        temp_output_dir: str = "/scratch"):
    """
    Given a list of files, process them one at a time. Gather the files written 
    during processing and return that list. 

    Parameters
    ----------
        file_path_list
            list of strs or pathlib.Path objs, assumed to be associated with
            gzipped EMBL/GenBank flat files. 
        database_params 
            dict or configparser.ConfigParser obj. necessary keys or 
            attributes are "user", "password", "host", and "port".
        db_name
            str, name of the EFI database to be used to perform queries.
        final_output_dir
            str, global path for storage space on the HPC filesystem within 
            which result files will be moved to.
        temp_output_dir
            str, global path for storage space on the compute resource within
            which result files will be written. A temporary space for fast IO.
            Default = "/scratch"

    Returns
    -------
        "process_many_files"
            str, used to ID type of task.
        final_tab_files
            list of strs, corresponding to the tsv files written during the
            task, using the final_output_dir path.
        `time.time() - st`
            float, elapsed time for this task, units: seconds.
        file_path_list
            str, same as given input.

    """
    st = time.time()
    # use regex to match the parent directories' names; three layers worth if
    # in `wgs` tree of ENA or two layers worth if in `sequence` tree. This
    # regex will match a file path string, creating a list of a tuple with len
    # 5. First three elements are associated with the wgs tree, the remaining
    # two with the sequence tree. 
    # NOTE: THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    dir_pattern = re.compile(r"(wgs)\/(\w*)\/(\w*)|(sequence)\/(\w*)")
    # use regex to match the file name stem from the given file path; will 
    # create a list of len 1. 
    file_pattern = re.compile(r"\/(\w*)\.dat\.gz")

    # apply the regex on the first file string in file_path_list, only grab 
    # groups that were successfully matched. 
    # NOTE: this assumes that all files in the file_path_list are sourced from
    # the same directory; this will be a bug if files from different source dirs
    # are included in file_path_list
    matches = [elem for elem in dir_pattern.findall(file_path_list[0])[0] if elem]
    # create an output_dir string that easily maps to the files being parsed.
    # format will be e.g. "wgs-public-wds" or "sequence-con"
    if temp_output_dir:
        if temp_output_dir[-1] != "/":
            temp_output_dir += "/"
        out_dir = temp_output_dir + "-".join(matches)
        # make the directory
        os.makedirs(out_dir, exist_ok=True)
    else:
        if final_output_dir[-1] != "/":
            final_output_dir += "/"
        out_dir = final_output_dir + "-".join(matches)
        # make the directory
        os.makedirs(out_dir, exist_ok=True)
   
    # connect to the database
    db_connection = mysql_database.IDMapper(database_params, db_name)

    tab_files = []
    for file_path in file_path_list:
        start_time = time.time()
        # grab the stem of the file name to use in writing results
        fn_name = file_pattern.findall(file_path)[0]
        tab_file = out_dir + f"/{fn_name}.tab"
        # process the file
        parse_embl.process_file(
            file_path, 
            db_connection, 
            tab_file
        )
        stop_time = time.time()
        # if the process_file() call does not output any results, no tab_file 
        # will be written so check to see if the expected file exists.
        if os.path.isfile(tab_file):
            tab_files.append(tab_file)
    
    db_connection.close()
    
    # if tab_files is empty, no files need to be shutil'd from temp to final
    # storage spaces.
    if not tab_files:
        final_tab_files = [""]
    # if temp_output_dir was actually used, need to shutil.move all the files
    # written to a scratch space to the final_output_dir location
    elif temp_output_dir:
        final_tab_files = []
        if final_output_dir[-1] != "/":
            final_output_dir += "/"
        # follow the same directory naming scheme as used in the temp space
        # NOTE: this assumes that all files in the file_path_list are sourced
        # from the same directory; this will be a bug if files from different
        # source dirs are included in file_path_list
        final_dir = final_output_dir + "-".join(matches)
        os.makedirs(final_dir, exist_ok=True)
        # loop over tab files and move them from the temp to the final storage
        # space; shutil.move() returns the new path string for the moved file
        for tab_file in tab_files:
            new_tab_file = shutil.move(tab_file, final_dir)
            final_tab_files.append(new_tab_file)
    # temp_output_dir was not used, so no shutil'ing needs to be done; tab 
    # files were already written to their final destination
    else:
        final_tab_files = tab_files
    
    return "process_many_files", final_tab_files, time.time() - st, file_path_list


