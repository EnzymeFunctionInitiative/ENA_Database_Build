
import re
import gzip

###############################################################################
# Parsing an EMBL flat file
###############################################################################

def process_file(
        file_path: str, 
        database_connection,
        output_file: str):
    """
    read gzipped GenBank or EMBL file, loop over entries, search for lines 
    specifying sequences and their EMBL IDs as well as Uniprot IDs and others 
    (if available).

    PARAMETERS
    ----------
        file_path
            string or pathlib.Path, assumed to be a gzipped embl flatfile.
        database_connection
            IDMapper object, connection object to the MYSQL database that 
            will be queried. 
        output_file
            string or pathlib.Path, file to be written with results from the
            function

    RETURN
    ------
        output_file
            string or pathlib.Path, file written to with results; if no 
            results were found, this file will not actually be written so a 
            check for its existence outside of this function is necessary
    """
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

    # create the dict to be used to gather all parsing results
    process_results = {}
    # instantiate vars to be filled; 
    ID = ""
    uniprotIds = []
    proteinIds = []
    count = 0
    CHR = 0
    DIR = 0
    START = 0
    END = 0

    # open and read the gzipped file
    with gzip.open(file_path, 'rt') as f:
        # loop over each line in f without reading the whole file
        for line in f: 
            # apply the search regex on line
            search_results = search_pattern.findall(line)
            
            # if the line does not match the regex pattern, move on to the 
            # next line
            if not search_results:
                continue
            
            # the line does match the regex pattern
            # turn search_results into the tuple of len 4
            search_results = search_results[0]
            
            # regex search found a new ENA Accession ID and DNA/chromosome
            # type
            if search_results[0] and search_results[1]:
                # need to handle the previous entry's data
                # check that the last count value hasn't been added to the
                # previous ID's subdir AND count != 0
                if count and START and END:
                    # the results for the file's last ID need to be added
                    # to the process_results dict
                    process_results[ID][count] = {
                        "uniprotIds": uniprotIds,
                        "proteinIds": proteinIds,
                        "seqcount": count,
                        "CHR": CHR,
                        "DIR": DIR,
                        "START": START,
                        "END": END,
                    }
                    # empty the lists and reassign values
                    uniprotIds = []
                    proteinIds = []
                    START = 0
                    END = 0
                    CHR = 0
                    DIR = 0
                    count = 0

                # now handle the new line's data
                # search_results[0] is the ENA Accession ID
                ID = search_results[0]
                # search_results[1] is the type of chromosome structure, 
                # linear or circular; there are some non-standard 
                # structures
                if search_results[1] in ["linear","circular"]:
                    # check if the type is either linear or circular
                    CHR = 1 if search_results[1] == "linear" else 0
                else:
                    print(f"!!! Unknown chromosome type observed in {file_path}: {line}")
                    # by replacing ID with an empty string, we are 
                    # effectively ignoring unexpected chromosome type 
                    # strings; we'll remove the "" key from the dict in
                    # the end
                    ID = ""
                
                # create the dict entry in process_results, key being the
                # ENA accession ID and value being a dict
                process_results[ID] = {}
                # this subdict will be filled with CDS entries, keys being
                # the count and values being dict of "uniprotIDs", 
                # "proteinIDs", "CHR", "DIR", "START", and "END"
                count = 0   # haven't found any CDS yet

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
                # END are non-zero, then we've found a new CDS line after
                # already having collected the previous CDS's results. 
                if count and START and END:
                    # create the associated subdict to the ID's dict
                    process_results[ID][count] = {
                        "uniprotIds": uniprotIds,
                        "proteinIds": proteinIds,
                        "seqcount": count,
                        "CHR": CHR,
                        "DIR": DIR,
                        "START": START,
                        "END": END,
                    }
                    # empty the lists
                    uniprotIds = []
                    proteinIds = []
                
                # found a CDS line denoting a coding sequence, add one to
                # the sequence count
                count += 1

                # determine the directionality of the CDS
                if "complement" in line:
                    DIR = 0
                else:
                    DIR = 1

                # find the start and stop values for the sequence
                cds_result = cds_pattern.findall(line)
                if cds_result:
                    START, END = cds_result[0]

    # check that the last count value hasn't been added to the 
    # process_results[ID] subdir AND count != 0
    if count and START and END:
        # the results for the file's last ID need to be added to 
        # the process_results dict
        process_results[ID][count] = {
            "uniprotIds": uniprotIds,
            "proteinIds": proteinIds,
            "seqcount": count,
            "CHR": CHR,
            "DIR": DIR,
            "START": START,
            "END": END,
        }

    # remove the false entry subdict if it exists
    if "" in process_results.keys():
        process_results.pop("")
   
    ## past version of code did not use processed_already
    #processed_already = {}

    # now do the reverseLookup against the SQL database
    # loop over each key in the process_results dict, loop over loci 
    # associated with that ena_id, and consider the proteinIds 
    # list as foreign_ids in a IDMapper.reverse_lookup call
    for ena_id in process_results.keys():
        for loci in process_results[ena_id].keys():
            loci_subdict = process_results[ena_id][loci]
            # perform the reverse lookup
            rev_uniprot_ids, no_match = database_connection.reverse_lookup(loci_subdict["proteinIds"])
            
            ## filter rev_uniprot_ids to only include Ids not already in processed_already
            #rev_uniprot_ids_to_add = [uniprot_id for uniprot_id in rev_uniprot_ids if uniprot_id not in processed_already]
            
            # check whether the rev_uniprot_ids set is empty
            if not rev_uniprot_ids:
                # if it is empty, use the 
                # process_results[ena_id][loci]["uniprotIds"] list; maybe 
                # should be the combined set of the two lists?
                uniprot_ids = loci_subdict["uniprotIds"]
            else:
                uniprot_ids = rev_uniprot_ids
            
            # loop over uniprot_ids and append to file. If uniprot_ids is
            # empty, no file is written.
            for id_ in uniprot_ids:
                with open(output_file,"a") as out_tab:
                    out_tab.write(f"{ena_id}\t{id_}\t{loci_subdict['seqcount']}\t{loci_subdict['CHR']}\t{loci_subdict['DIR']}\t{loci_subdict['START']}\t{loci_subdict['END']}\n")

    return output_file


