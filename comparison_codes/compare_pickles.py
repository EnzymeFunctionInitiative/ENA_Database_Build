
import sys
import pickle
import logging

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
# Comparison function
###############################################################################

def compare_two_dicts_of_dicts(dict1, dict1_name, dict2, dict2_name, log_file):
    """
    code to loop over all keys, subkeys, and subsublists in dict1 to see if 
    they are present in dict2. If not, log the absence/mismatch between the two
    dicts. 

    Parameters
    ----------
        dict1,
            dict of dict of lists. 
        dict1_name,
            str, textual description of dict1 used in logging.
        dict2,
            dict of dict of lists. 
        dict2_name,
            str, textual description of dict2 used in logging.
        log_file,
            str or pathlib.Path, path to where the log file will be created

    Return
    ------
        missing_keys,
            int, number of keys present in dict1 but absent in dict2
        missing_subkeys,
            int, number of subkeys present in dict1's subdictionarys that are 
            absent in dict2's subdictionaries
        mismatches,
            int, number of subsublists that do not exactly match between dict1 
            and dict2
        mismatches_list,
            list of ints, len = 5, number of mismatch instances that map to the 
            subsublist element(s)

    """
    # create the variables to be filled
    missing_keys = 0
    missing_subkeys  = 0
    mismatches = 0
    mismatches_list = [0,0,0,0,0]
    nComparisons = 5
   
    # start the logger object
    dict_comp_logger = setup_logger("dict_comparison", log_file)
    # loop over dict1 keys, get the associated subdict and then check if the
    # key:subdict pair exists in dict2
    for key in dict1.keys():
        dict1_subdict = dict1[key]
        dict2_subdict = dict2.get(key, {})
        # if the key doesn't map to a subdict, then log it and move on
        if not dict2_subdict:
            missing_keys += 1
            dict_comp_logger.info(f"{key} missing from {dict2_name}, but in {dict1_name}.")
            continue
        
        # loop over subkeys in the dict1_subdict, get the associated subsublist
        # and check if the subkey:subsublist pair exists in dict2's subdict
        for subkey in dict1_subdict.keys():
            dict1_subsublist = dict1_subdict[subkey]
            dict2_subsublist = dict2_subdict.get(subkey, [])
            # if the subkey doesn't map to a subsublist, then log it and move on
            if not dict2_subsublist:
                missing_subkeys += 1
                dict_comp_logger.info(f"{key}, {subkey} missing from {dict2_name}, but in {dict1_name}.")
                continue
            
            # both key and subkey are in dict1's and dict2's archetecture. Now 
            # check that the subsublist elements match. Since we're interested
            # in knowing which element(s) did not match, compare each element 
            # individually.
            matches = [dict1_subsublist[i] == dict2_subsublist[i] for i in range(nComparisons)]
            # if there is a mismatch, determine which element and log it
            if sum(matches) != nComparisons:
                #dict_comp_logger.info(f"{key}, {subkey}, {dict1_name}: {dict1_subsublist}, {dict2_name}: {dict2_subsublist}")
                mismatches += 1
                for i in range(nComparisons):
                    mismatches_list[i] += 1 if not matches[i] else 0
                    if i > 0:
                        dict_comp_logger.info(f"{key}, {subkey}, {dict1_name}: {dict1_subsublist}, {dict2_name}: {dict2_subsublist}")
    
    # close the logger object
    clean_logger(dict_comp_logger)
    
    return missing_keys, missing_subkeys, mismatches, mismatches_list


###############################################################################
# Main
###############################################################################

if __name__ == "__main__":
    # assign input args to variables
    perl_pkl = sys.argv[1]
    pyth_pkl = sys.argv[2]
    log_dir  = sys.argv[3]

    # start up the main logger object
    main_logger = setup_logger("comparison_logger",log_dir + "/main.log")
    main_logger.info("Loading up the pickle files.")
    
    # read the first pickle file to gather the stored dict of dict of lists
    with open(perl_pkl,'rb') as perl_tab:
        contents, perl_dict = pickle.load(perl_tab)
        main_logger.info(f"Done loading the {perl_pkl} file.")
    
    # read the second pickle file to gather the stored dict of dict of lists
    with open(pyth_pkl,'rb') as pyth_tab:
        contents, pyth_dict = pickle.load(pyth_tab)
        main_logger.info(f"Done loading the {pyth_pkl} file.")
    
    # log and run the comparison
    main_logger.info('#'*80)
    main_logger.info('Starting comparison of the Perl ENA tab to the Python ENA tab.')
    python_vs_perl_logger = log_dir + "/python_vs_perl.log"
    missing_protIds, missing_uniIds, mismatches, mismatches_list = compare_two_dicts_of_dicts(pyth_dict, "Python ENA tab",perl_dict, "Perl ENA tab", python_vs_perl_logger)
    main_logger.info(f'Perl ENA tab is missing {missing_protIds:,} Record entries, relative to the Python ENA tab.')
    main_logger.info(f'Perl ENA tab is missing {missing_uniIds:,} locus entries, relative to the Python ENA tab.')
    main_logger.info(f'There are {mismatches:,} mismatches between locus entries in the two files. Specifically:\n\t{mismatches_list[0]:,} disagreements in the seqCount column\n\t{mismatches_list[1]:,} disagreements in the CHR column\n\t{mismatches_list[2]:,} disagreements in the DIR column\n\t{mismatches_list[3]:,} disagreements in the START column\n\t{mismatches_list[4]:,} disagreements in the END column')

    # log and run the reverse comparison
    main_logger.info('#'*80)
    main_logger.info('Starting comparison of the Python ENA tab to the Perl ENA tab.')
    perl_vs_python_logger = log_dir + "/perl_vs_python.log"
    missing_protIds, missing_uniIds, mismatches, mismatches_list = compare_two_dicts_of_dicts(perl_dict, "Perl ENA tab",pyth_dict, "Python ENA tab", perl_vs_python_logger)
    main_logger.info(f'Python ENA tab is missing {missing_protIds:,} Record entries, relative to the Perl ENA tab.')
    main_logger.info(f'Python ENA tab is missing {missing_uniIds:,} locus entries, relative to the Perl tab.')
    main_logger.info(f'Number of mismatches already reported above.')
    main_logger.info('#'*80)
    
    # close the logger object
    clean_logger(main_logger)


