
import time
import glob
import gzip
import re
import os
import shutil

import mysql_database
import parse_embl
import bio_parse_embl

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
    # to us. As far as I know, the second underscored section of the file name
    # denote the origin species type, which is what we need to consider.
    # NOTE: THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    if "sequence" in dir_path:
        # NOTE: regex to only gather file names with (ENV|PRO|FUN|PHG) in them
        pattern = re.compile(r"_(ENV|PRO|FUN|PHG)_")
        files = [file_ for file_ in files if pattern.search(file_)]

    return "glob_files", files, time.time() - st, dir_path
    ## could gather file size as well as name to enable sorting the files from
    ## largest to smallest; potential to optimize task prioritization
    #return "glob_files", [(file, os.path.getsize(file)) for file in files], time.time() - st, dir_path


def process_many_files(
        file_path_list: list, 
        database_params: dict,
        db_name: str, 
        final_output_dir: str,
        temp_output_dir: str = "/scratch"):
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
        final_output_dir
            string, global path for storage space on the HPC filesystem within 
            which result files will be moved to.
        temp_output_dir
            string, global path for storage space on the compute resource within
            which result files will be  written. A temporary space for fast IO.
            Default = "/scratch"

    RETURNS
    -------
        "process_many_files"
            string used to ID type of task
        final_tab_files
            list of strings corresponding to the tsv files written during the
            task, using the final_output_dir path
        `time.time() - st`
            elapsed time for this task, units: seconds
        file_path_list
            same as given input

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
    file_pattern= re.compile(r"\/(\w*)\.dat\.gz")

    # apply the regex on the first file string in file_path_list, only grab 
    # groups that were successfully matched. 
    # NOTE: this assumes that all files in the file_path_list are sourced from
    # the same directory; this will be a bug if files from different source dirs
    # are included in file_path_list
    matches = [elem for elem in dir_pattern.findall(file_path_list[0])[0] if elem]
    # create an output_dir string that easily maps to the files being parsed
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
        out_diir = final_output_dir + "-".join(matches)
        # make the directory
        os.makedirs(out_dir, exist_ok=True)
   
    # connect to the database
    db_connection = mysql_database.IDMapper(database_params, db_name)

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
        # if the file does not return any results, no file will be written so 
        # check to see if the expected file exists.
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


