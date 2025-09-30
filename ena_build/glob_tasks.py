
import re
import os
from typing import List, Tuple

###############################################################################
# Functions used as Dask Tasks
###############################################################################

def glob_subdirs(dir_path: str) -> Tuple[str, List[str], float, str]:
    """
    Search for subdirectories in the provided directory path.

    Parameters
    ----------
    dir_path: str 
        global or local path within which the search for subdirs will occur.

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


def glob_files(
        dir_path: str,
        file_filter_pattern: None | re.Pattern = None
    ) -> Tuple[str, List[str], float, str]:
    """
    Return list of files matching the search string.
    
    Parameters
    ----------
        dir_path
            str, global or local path within which the search for subdirs will
            occur.
        file_filter_pattern
            None or re.Pattern object, a text pattern used to filter files from
            the list of files produced by os.scandir(). Defualt: None

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
   
    # apply the file filter pattern
    if file_filter_pattern:
        # filter files based on whether they match the file_filter_pattern
        files = [_file for _file in files if file_filter_pattern.search(_file)]

    return "glob_files", files, time.time() - st, dir_path



    ## Only a subset of data files in the ENA sequence/ subdir are of interest 
    ## to us. As far as I know, the second underscored section of the file name
    ## denote the origin species type, which is what we need to consider.
    ## NOTE: THIS MAY BE A BUG DEPENDING ON CHANGES MADE BTW ENA VERSIONS
    #if "sequence" in dir_path:
    #    # NOTE: regex to only gather file names with (ENV|PRO|FUN|PHG) in them
    #    pattern = re.compile(r"_(ENV|PRO|FUN|PHG)_")
    #    files = [file_ for file_ in files if pattern.search(file_)]

