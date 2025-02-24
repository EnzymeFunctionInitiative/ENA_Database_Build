
import sys
import pickle
import logging

perl_pkl = sys.argv[1]
pyth_pkl = sys.argv[2]
log_dir  = sys.argv[3]

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
# Main
###############################################################################
if __name__ == "__main__":

    main_logger = setup_logger("comparison_logger",log_dir + "/main.log")
    main_logger.info("Loading up the pickle files.")
    
    with open(perl_pkl,'rb') as perl_tab:
        contents, perl_dict = pickle.load(perl_tab)
        main_logger.info(f"Done loading the {perl_pkl} file.")
    
    with open(pyth_pkl,'rb') as pyth_tab:
        contents, pyth_dict = pickle.load(pyth_tab)
        main_logger.info(f"Done loading the {pyth_pkl} file.")
    
    main_logger.info('#'*80)
    main_logger.info('Starting comparison of the Perl tab file to the Python tab file.')
    missing_protIds = 0
    missing_uniIds  = 0
    mismatches = 0

    python_vs_perl_logger = setup_logger("py_vs_prl", log_dir + "/python_vs_perl.log")
    # loop over pyth_dict keys, get the subdicts, and compare each key in the 
    # gathered subdicts
    for proteinId in pyth_dict.keys():
        pyth_subdict = pyth_dict[proteinId]
        perl_subdict = perl_dict.get(proteinId, {})
        if not perl_subdict: 
            missing_protIds += 1
            python_vs_perl_logger.info(f"{proteinId} missing from the Perl ENA tab, but in the Python ENA tab")
            continue
        
        #UniProtIds = set()
        # loop over subdict keys in both subdicts
        for uniprot in pyth_subdict.keys():
            ## hopefully there aren't instances of duplicate UniProtIds in the chromosome's subdict
            #if uniprot in UniProtIds:
            #    continue
            #UniProtIds.add(uniprot)
            # gather the two lists
            pyth_subsublist = pyth_subdict[uniprot]
            perl_subsublist = perl_subdict.get(uniprot,[])
            
            # check whether perl_tab had the expected list
            if not perl_subsublist:
                missing_uniIds  += 1
                python_vs_perl_logger.info(f"{proteinId}, {uniprot} missing from the Perl ENA tab, but in the Python ENA tab.")
                continue
            
            # count the number of times where the two subsublists match
            matches = sum([pyth_subsublist[i] == perl_subsublist[i] for i in range(len(pyth_subsublist))])
            # if any of the elements don't match, print about it
            if matches != len(pyth_subsublist):
                #python_vs_perl_logger.info(f"{proteinId}, {uniprot}, python: {pyth_subsublist}, perl: {perl_subsublist}")
                mismatches += 1
    
    main_logger.info(f'Perl tab file is missing {missing_protIds:,} Record entries, relative to the Python tab file.')
    main_logger.info(f'Perl tab file is missing {missing_uniIds:,} locus entries, relative to the Python tab file.')
    main_logger.info(f'Python tab file has {mismatches:,} mismatches between locus entries with the Perl tab file.')
    
    main_logger.info('#'*80)
    main_logger.info('Starting comparison of the Python tab file to the Perl tab file.')
    missing_protIds = 0
    missing_uniIds  = 0
    mismatches = 0
    
    perl_vs_python_logger = setup_logger("prl_vs_py", log_dir + "/perl_vs_python.log")
    # loop over perl_dict keys, get the subdicts, and compare each key in the 
    # gathered subdicts
    for proteinId in perl_dict.keys():
        perl_subdict = perl_dict[proteinId]
        pyth_subdict = pyth_dict.get(proteinId, {})
        if not pyth_subdict: 
            missing_protIds += 1
            perl_vs_python_logger.info(f"{proteinId} missing from the Python ENA tab, but in the Perl ENA tab")
            continue
        
        #UniProtIds = set()
        # loop over subdict keys in both subdicts
        for uniprot in perl_subdict.keys():
            ## hopefully there aren't instances of duplicate UniProtIds in the chromosome's subdict
            #if uniprot in UniProtIds:
            #    continue
            #UniProtIds.add(uniprot)
            # gather the two lists
            perl_subsublist = perl_subdict[uniprot]
            pyth_subsublist = pyth_subdict.get(uniprot,[])
            
            # check whether pyth_tab had the expected list
            if not pyth_subsublist:
                missing_uniIds  += 1
                perl_vs_python_logger.info(f"{proteinId}, {uniprot} missing from the Python ENA tab, but in the Perl ENA tab.")
                continue
            
            # count the number of times where the two subsublists match
            matches = sum([pyth_subsublist[i] == perl_subsublist[i] for i in range(len(perl_subsublist))])
            # if any of the elements don't match, print about it
            if matches != len(perl_subsublist):
                #perl_vs_python_logger.info(f"{proteinId}, {uniprot}, python: {pyth_subsublist}, perl: {perl_subsublist}")
                mismatches += 1

    main_logger.info(f'Python tab file is missing {missing_protIds:,} Record entries, relative to the Perl tab file.')
    main_logger.info(f'Python tab file is missing {missing_uniIds:,} locus entries, relative to the Perl tab.')
    main_logger.info(f'Perl tab file has {mismatches:,} mismatches between locus entries with the Python tab file.')
    
    main_logger.info('#'*80)
    
    clean_logger(main_logger)
    clean_logger(perl_vs_python_logger)
    clean_logger(perl_vs_python_logger)


