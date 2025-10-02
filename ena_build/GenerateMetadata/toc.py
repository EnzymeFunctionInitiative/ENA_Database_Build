
import os
import hashlib
import gzip

###############################################################################
# Functions to gather metadata about a file
###############################################################################

def md5_of_gzip_file(file_path: str) -> str:
    """
    Calculates the MD5 hash of a gzip file.

    Arguments
    ---------
        file_path
            str, the path to the gzip file.

    Returns
    -------
        The MD5 hash of the file as a hexadecimal string.
    """
    md5_hash = hashlib.md5()
    with gzip.open(file_path, 'rb') as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def get_file_stats(
        file_path: str,
        key_list: List[str]
    ) -> Dict[str, Any]:
    """
    Gets the stats associated with the file at file_path.

    Arguments
    ---------
        file_path
            str, the path to the file.
        key_list
            list of str, os.stat_result attribute names to be gathered from
            the os.stat() call on file_path. Defaults to ["st_mtime"] (last
            modified time; reported in epoch seconds).

    Returns
    -------
        dict, relevant metadata information about the file. Which keys is
        dependent on the key_list input argument.
    """
    stat_result = os.stat(file_path)
    return {key: stat_result.__getattribute__(key) for key in key_list}


###############################################################################
# Wrapper function for the above funcs
###############################################################################

def get_metadata(
        file_path: str,
        md5_bool: bool = False,
        key_stats_list: List[str] = ["st_mtime"]
    ) -> Dict[str,Any]:
    """
    Wrapper function for the get_file_stats() and md5_of_gzip_file() functions.

    Arguments
    ---------
        file_path
            str, the path to the file.
        md5_bool
            bool, controls whether the md5 hash is calculated for the given
            file.
        key_stats_list
            list of str, strings must match the attribute names of a
            `os.stat_result` object. See https://docs.python.org/3/library/os.html#os.stat_result
            for documentation.

    Returns
    -------
        dict, contains the metadata values associated with the key_stats_list
        attributes as well as a md5_hash key:value pair (if md5_bool == True).
    """
    # NOTE: this list should originate from a TOC sql table's columns and should
    # be passed from the workflow code to the task function to this function
    # instead of being hard-defined here.
    key_stat_list = [
        "st_size",  # units of bytes
        "st_mtime" # units of epoch seconds
    ]

    # verify the key list being used to gather metadata.
    if type(key_stats_list) != list:
        raise ValueError("User specified os.stat() keys to be gathered as"
            + " metadata is incorrectly formatted")

    file_stats = get_file_stats(file_path, key_stats_list)

    if md5_bool:
        file_stats.update(
            {
                "md5_hash": md5_of_gzip_file(file_path)
            }
        )

    return file_stats

