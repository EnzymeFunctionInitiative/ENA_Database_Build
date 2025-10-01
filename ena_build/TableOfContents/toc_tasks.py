
import time
import os
from typing import List, Tuple, Dict, Any
from dataclass import dataclass, field

import toc
import mapping

###############################################################################
# Data Class for Gathering File Metadata
###############################################################################

@dataclass(repr = False, eq = False, match_args = False)
class FileMetadata:
    """
    Class for keeping track of metadata associated with a specific file.

    Attributes
    ----------
        file_path
            str, path to the file associated with this object instance.
        total_processing_time
            float, units: seconds. Amount of time spent processing the file.
        toc
            dict, table of contents metadata dictionary. Contents of this dict
            depend on the toc.get_metadata() function. May be an empty dict.
        ids
            list, protein_ids found within the file. May be an empty list.

    Both `toc` and `ids` attributes can only be defined via calling their
    explicit keyword when a FileMetadata object is instantiated. Neither of
    these attributes are printed when __repr__() is called for this class.
    """
    file_path: str
    total_processing_time: float
    toc: Dict[str, Any] = field(kw_only = True, default_factor=dict)
    ids: List[str] = field(kw_only = True, default_factor=list)

    def __repr__(self):
        # keep the printed representation of the object simple since toc and
        # ids can be complex/large.
        return (f"FileMetadata(file_path={self.file_path},"
                + f" total_processing_time={self.total_processing_time})")


###############################################################################
# Functions used as Dask Tasks
###############################################################################

def gather_files_metadata(
        file_path_list: List[str],
        toc_bool: bool = False,
        md5_hash_bool: bool = False,
        mapping_bool: bool = False,
    ) -> Tuple[str, List[FileMetadata], float]:
    """
    Given a list of files, process them one at a time. Gather metadata (if
    toc_bool is True) and/or protein_ids found within the file (if mapping_bool
    is True). Return a list of FileMetadata class objects.

    Parameters
    ----------
        file_path_list
            list of strs or pathlib.Path objs, assumed to be associated with
            gzipped EMBL/GenBank flat files.
        toc_bool
            bool, if True, gather metadata about the files in file_path_list.
            Default: False.
        md5_hash_bool
            bool, if True, determine the md5sum hash for the files in
            file_path_list. Default: False.

    Returns
    -------
        "gather_files_metadata"
            str, used to ID type of task.
        metadata_obj_list
            list, list of FileMetadata object instances that contain the
            important metadata for each file in file_path_list.
        `time.time() - st`
            float, elapsed time for this task, units: seconds.

    """
    st = time.time()
    metadata_obj_list = []
    for file_path in file_path_list:
        start_time = time.time()

        # gather TOC information
        if toc_bool:
            toc_contents = toc.get_metadata(file_path, md5_hash_bool)
        # or not
        else:
            toc_contents = {}

        # gather protein_id mapping information
        if mapping_bool:
            ids_list = mapping.process_file(file_path)
        # or not
        else:
            ids_list = []

        # stash the file's metadata (toc and/or ids) in the data class object
        # then stash that object instance in the list to be returned by the
        # function
        metadata_obj_list.append(
            FileMetadata(
                file_path,
                time.time() - start_time,
                toc = toc_contents,
                ids = ids_list
            )
        )

    return "gather_files_metadata", metadata_obj_list, time.time() - st


