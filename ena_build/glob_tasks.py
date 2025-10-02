
import re
import os
from typing import List, Tuple

###############################################################################
# Define regex pattern constants variables
###############################################################################

# SOURCE_PATTERN is used as a file filter to only consider files from the given
# sources
SOURCE_PATTERN = re.compile(r"_(ENV|PRO|FUN|PHG)_")

# DIR_PATTERN is used to parse subdirectories' names; three layers worth if in
# `wgs` tree of ENA or two layers worth if in `sequence` tree. When called via
# re.findall(), this regex will creat a list (len = 1) with a tuple with len 3.
# NOTE: THIS IS HIGHLY DEPENDENT ON THE DIRECTORY TREE STRUCTURE OF THE ENA
# DOWNLOAD
DIR_PATTERN = re.compile(r"(wgs|sequence)\/(\S*)\/(\S*)\/")
# NOTE: \S is a very greedy regex pattern; we should avoid using it...

# FILE_NAME_PATTERN is used to get the file name stem from the given path; will
# create a list of len 1.
# NOTE: this assumes that the stem of dat.gz files of interest only contain
# alphanumeric characters and underscores.
FILE_NAME_PATTERN = re.compile(r"\/(\w*)\.dat\.gz")

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


