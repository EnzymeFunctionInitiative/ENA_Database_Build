
###############################################################################
# CODE OVERVIEW of `buildena.nf`:
#
# ENA database is organized as follows: 
# > module load ena
# "ena/20241205 database is located at /private_stores/mirror/ena/20241205"
# > ls /private_stores/mirror/ena/20241205
#   sequence/ wgs/
#
# The buildena.nf code specifically looks for:
# /private_stores/mirror/ena/20241205/wgs/public/**/*dat.gz
# /private_stores/mirror/ena/20241205/wgs/suppressed/**/*dat.gz
# /private_stores/mirror/ena/20241205/sequence/**/*dat.gz # filters sequence files by { it =~ /(ENV|PRO|FUN|PHG)/ }
#       using regular expression: =~, True if given pattern occurs anywhere in the string
#       given pattern is any of ENV, PRO, FUN, PHG
#
# Files matching the glob and regex are concatenated into a list, 
#    ordered ~/wgs/public/**/*dat.gz, ~/sequence/**/*dat.gz (regex applied), ~/wgs/suppressed/**/*dat.gz
#
# File names are used to create tuples for each file, containing
#       file stem (removing .dat.gz from the end),
#       path to the file (including the .dat.gz),
#       f.parent.parent.baseName + "/" + f.parent.baseName... e.g. "public/aaa" from /private_stores/mirror/ena/20241205/wgs/public/abz/ABZX01.dat.gz
#
#
# The files were individually decompressed via gunzip and saved to storage (baaaaddddd) in a nextflow process (baaaaadddddd)
#
# The decompressed files were then parsed individually (baaaaadddddd) in a nextflow process (baaaaddddd)
# input to the parser is: 
#       filename, string, used for output file name
#       dat_file, path, points to the input data file path
#       dir_name, string to append to the end of ${params.output_dir}/
#
# parse_embl() process is run on each file's tuple... more info needed as to what this is grabbing from the unzipped files.
#
# ideas for the dask or nextflow pipeline: 
#   spin up three instances of a process... Channel.fromPath processes for the three directory paths, getting the list of subdirectories
#   given lists of subdirectories, spin up a process for each that aggregates the *.dat.gz files in the respective subdir
#   given the lists of *dat.gz files, shard the list of files into N even sets of files
#   spin up N processes that parse and output the data for each file in a given set; 
#       do not unzip to GPFS storage; if an gunzip call is absolutely necessary, use on node scratch or /dev/shm space to do so
#
###############################################################################

import os
import time
import glob
import gzip
import re
import argparse
import logging
import configparser
import uuid

#import platform
#import io
from pathlib import Path

#import dask
#import dask.config
from distributed import Client, Worker, as_completed, get_worker

import mysql.connector

from Bio import GenBank

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
        fn_name = Path(file_).name.split('.')[0]
        tab_file = process_file(file_, db_connection, out_dir + f"/{fn_name}.tab")
        stop_time = time.time()
        # if the file_ does not return tab results, no file will be written
        if os.path.isfile(out_dir + f"/{fn_name}.tab"):
            tab_files.append(tab_file)
    
    db_connection.close()
    
    # quick and dirty handling of the task output
    if not tab_files:
        tab_files = [""]

    return "process_many_files", tab_files, time.time() - st, file_path_list


###############################################################################
# Functions
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


def process_file(
        file_path: str, 
        database_connection,
        output_file: str):
    """
    read gzipped GenBank or EMBL file, loop over entries, search for lines 
    specifying sequences and their EMBL IDs as well as their Uniprot IDs (if 
    available).
    """
    # create the dict to be used to gather all parsing results
    process_results = {}
    # initialize vars used to gather info
    ID  = ""
    uniprotIds = []
    proteinIds = []
    count = 0
    CHR = ""
    DIR = -9999
    START = 0
    END   = 0

    # create the regex search strings before hand to improve efficiency of doing
    # the regex searches on each line of the file(s).
    search_strs = [
        r"(?:^ID\s+(\w+);\s\w\w\s\w;\s(\w+);\s.*)", # search for ID lines, group 0 and 1 map to ID and type of genome (circular or linear);
        r'(?:^FT\s+\/protein_id=\"([a-zA-Z0-9\.]+)\")', # search for "protein_id" FT lines, group 2 maps to the quoted ID. 
        r'(?:^FT\s+\/db_xref=\"UniProtKB\/[a-zA-Z0-9-]+:(\w+)\")',  # search for UniProtKB/... database accession ID lines, group 3 matches the associated accession ID. 
        r"(?:^FT\s+CDS)"    # search for "FT   CDS" lines, need to do boolean checks for the various lines following this pattern
    ]
    
    # compile the combined search pattern.
    # for each line that successfully matches one of the search strings, a list
    # of one tuple of len 4 will be created. The zeroth element of the tuple
    # maps to the ENA ID. First maps to the type of genome. Second maps to a
    # protein_id string. And third maps to a UniProtKB accession ID. Depending 
    # on which line is matched, some or all of these elements could be empty
    # strings _but_ the tuple will always be len 4. If the line does not match
    # any search strings, then an empty list will be returned. 
    search_pattern = re.compile("|".join(search_strs))
   
    # create the "FT   CDS" search pattern. 
    # If the line matches this search string, then the a list of a tuple of 
    # len 2 will be returned. The zeroth and first elements will be the START
    # and END values for the sequence, respectively. Else, the regex search 
    # will return an empty list. 
    cds_pattern = re.compile(r"(\d+)\..*\.\>?(\d+)")

    # open and read the gzipped file
    with gzip.open(file_path, 'rt') as f:
        # loop over each line in f without reading the whole file
        for line in f: 
            # apply the search regex on line
            search_results = search_pattern.findall(line)
            # if the line matches the regex pattern, parse it some more, else
            # move on to the next line
            if search_results:
                # turn search_results into the tuple of len 4
                search_results = search_results[0]
                
                # regex search found a new ENA Accession ID and DNA/chromosome
                # type
                if search_results[0] and search_results[1]:
                    
                    # before the new ID can be processed, add the previous 
                    # ID's results to the results dict.
                    # this will create one false entry at the start of the file
                    process_results[ID] = {
                        "uniprotIds": uniprotIds,
                        "proteinIds": proteinIds,
                        "seqcount": count,
                        "CHR": CHR,
                        "DIR": DIR,
                        "START": START,
                        "END": END,
                    }

                    # search_results[0] is the ENA Accession ID
                    ID = search_results[0]
                    # search_results[1] is the type of chromosome structure, 
                    # linear or circular
                    if search_results[1] in ["linear","circular"]:
                        # check if the type is either linear or circular
                        CHR = 1 if search_results[1] == "linear" else 0
                    #elif search_results[1] == "XXX":
                    #    print(f"!!! figure out what to do for {file_path}: {line}")
                    else:
                        print(f"!!! Unknown chromosome type observed in {file_path}: {line}")
                        # by replacing ID with empty string, we are effectively
                        # ignoring unexpected chromosome type strings; we'll 
                        # remove the "" key from the dict in the end
                        ID = ""
                    
                    # re-initialize variables as empty to be filled with 
                    # information
                    uniprotIds = []
                    proteinIds = []
                    count = 0
                    DIR = -9999
                    START = 0
                    END = 0
           
                # regex search found a protein_id line
                elif search_results[2]:
                    proteinIds.append(search_results[2])

                # regex search found a UniProtKB line
                elif search_results[3]:
                    uniprotIds.append(search_results[3])
                
                # regex search matched but tuple is filled with empties, must 
                # have found a "FT   CDS" line
                else:
                    # if the list of uniportIds is occupied and ints START and 
                    # END are non-zero, then we're done processing CDS lines. 
                    # We can skip to the next line and keep doing so until we 
                    # hit a new "ID" line.
                    if uniprotIds and START and END:
                        continue
                    
                    # add one to the sequence count
                    # should this be before the above boolean check? 
                    count += 1

                    # determine the directionality of the encoding sequence
                    # can "FT   CDS" lines have different directionality? 
                    # should we check to see if DIR has already been assigned a 
                    # value?
                    if "complement" in line:
                        DIR = 0
                    else:
                        DIR = 1

                    # check for start and stop values for the sequence
                    cds_matches = cds_pattern.findall(line)
                    if cds_matches:
                        # as above, should we check to see if START and END 
                        # have been assigned a non-zero value already?
                        START, END = cds_matches[0]

    # before moving on, the results for the file's last ID need to be added to 
    # the process_results dict
    process_results[ID] = {
        "uniprotIds": uniprotIds,
        "proteinIds": proteinIds,
        "seqcount": count,
        "CHR": CHR,
        "DIR": DIR,
        "START": START,
        "END": END,
    }

    # remove that first false entry before finishing
    process_results.pop("")
   
    ## past version of code did not use processed_already
    #processed_already = {}

    # now do the reverseLookup against the SQL database
    # loop over each key in the process_results dict and consider the proteinIds 
    # list as foreign_ids in a IDMapper.reverse_lookup call
    for ena_id in process_results.keys():
        # perform the reverse lookup
        rev_uniprot_ids, no_match = database_connection.reverse_lookup(process_results[ena_id]["proteinIds"])
        
        ## filter rev_uniprot_ids to only include Ids not already in processed_already
        #rev_uniprot_ids_to_add = [uniprot_id for uniprot_id in rev_uniprot_ids if uniprot_id not in processed_already]
        
        # check whether the rev_uniprot_ids list is empty
        if not rev_uniprot_ids:
            # if it is empty, use the process_results[ena_id]["uniprotIds"] 
            # list; maybe should be the combined set of the two lists?
            uniprot_ids = process_results[ena_id]["uniprotIds"]
        else:
            uniprot_ids = rev_uniprot_ids
        
        # loop over the list of uniprot_ids and append to file. If uniprot_ids
        # list is empty, no file is written.
        for id_ in uniprot_ids:
            with open(output_file,"a") as out_tab:
                out_tab.write(f"{ena_id}\t{id_}\t{process_results[ena_id]['seqcount']}\t{process_results[ena_id]['CHR']}\t{process_results[ena_id]['DIR']}\t{process_results[ena_id]['START']}\t{process_results[ena_id]['END']}\n")

    return output_file


class IDMapper:
    def __init__(self, database_params, db_name):
        """ Set up a connection to the user-specified MySQL database """
        try:
            cnx = mysql.connector.connect(
                user=database_params["user"],
                password=database_params["password"],
                host=database_params["host"],
                port=database_params["port"],
                database=db_name,
            )
        except mysql.connector.Error as err:
            if err.errno == mysql.connector.errorcode.ER_ACCESS_DENIED_ERROR:
                print(f"Something is wrong with your user name or password.\n{err}")
            elif err.errno == mysql.connector.errorcode.ER_BAD_DB_ERROR:
                print("Database {database_params.db_name} does not exist.\n{err}")
            else:
                print(err)
            cnx = ""
        
        self.dbh = cnx

    def close(self):
        """ Close the database connection """
        if self.dbh.is_connected():
            self.dbh.close()

    def reverse_lookup(self, foreign_ids, batch_size=1000):
        """
        """
        # make sure the foreign_ids list is not empty
        if not foreign_ids:
            return [], []
        
        # prep the sql string for the list of foreign_ids
        # this search is grabbing both foreign_id and uniprot_id from the table
        # if the row's foreign_id is in the list of foreign_ids (input).
        placeholder = ", ".join(["%s"]*len(foreign_ids))
        sql = f"SELECT foreign_id, uniprot_id FROM idmapping WHERE foreign_id IN ({placeholder})"
        
        # spin up a cursor on the database handle
        cursor = self.dbh.cursor(dictionary=True)
        # execute query: 
        # for each foreign_id found by the IN statement, this query will
        # produce a dict with keys "foreign_id" and "uniprot_id" w/ associated 
        # values. 
        # fetchone() will return this dict, fetchall() or fetchmany() will 
        # return a list of dicts associated with the successful queries
        # if a foreign_id is not found in the table but is in the list of 
        # foreign_ids, no dict will be returned. 
        cursor.execute(sql, foreign_ids)
        
        # prep the output arrays. will add found uniprot_ids to that list and
        # will remove the associated foreign_id from the no_matches list.
        uniprot_ids = []
        no_matches = set(foreign_ids)

        # instead of fetching results one or all at a time, let's grab a batch
        # with controlled size. This prevents inefficiencies of one at a time
        # while avoiding out of memory concerns if the query returns a huge 
        # number of results
        while True:
            # fetch a set of query results, number of rows dependent on
            # batch_size
            batch = cursor.fetchmany(batch_size)
            
            # if batch is an empty list, we've hit the end of query results so
            # break out of the while loop
            if not batch:
                break
            
            # loop over rows in the query results
            for row in batch:
                # add the uniprot_id to the list
                uniprot_ids.append(row['uniprot_id'])
                # remove the foreign_id from the no_matches list
                no_matches.discard(row['foreign_id'])
        
        # end cursor instances to keep dbh connections minimal
        cursor.close()
        
        return uniprot_ids, no_matches


