
import re
import gzip

###############################################################################
# Parsing an EMBL flat file
###############################################################################

class Record():
    def __init__(self, ID: str, CHR: int):
        """
        """
        self.ENA_ID = ID
        self.CHR: CHR,
        self.count: 0,
        self.uniprotIds: set(),
        self.proteinIds: set(),
        # this subdict will be filled with CDS entries, keys being the count 
        # and values being dict of "uniprotIDs", "proteinIDs", "DIR", "START",
        # and "END"
        self.loci_dict = {}
        self.current_locus = {
            "uniprotIds": set(),
            "proteinIds": set(),
            "DIR": -9999,
            "START": 0,
            "END": 0,
        }

    def add_locus(locus_id):
        """
        """
        self.loci_dict[locus_id] = self.current_locus
        self.current_locus = {
            "uniprotIds": set(),
            "proteinIds": set(),
            "DIR": -9999,
            "START": 0,
            "END": 0,
        }

    def process_record(last_locus, db_cnx, output_file):
        """
        """
        # check to see if the last_locus is not a key in the Record's loci_dict
        # as well as that locus' START, END, and proteinIds set values are 
        # equivalent to True
        if (last_locus not in self.loci_dict.keys()
                and self.current_locus["START"]
                and self.current_locus["END"]
                and self.current_locus["proteinIds"]):
            self.add_locus(last_locus)
        
        # loop over each locus in the Record object, doing a reverseLookup on
        # the locus' proteinIds to gather any uniprotIds. Write to file if the
        # locus has an associated uniprotId.
        for locus in self.loci_dict.keys():
            locus_subdict = self.loci_dict[locus]
            # perform the reverse lookup
            rev_uniprot_ids, no_match = db_cnx.reverse_lookup(
                list(locus_subdict['proteinIds'])
            )
            # check whether the rev_uniprot_ids set is empty
            if not rev_uniprot_ids:
                # if it is empty, use the loci's uniprotIds value
                # instead
                # maybe should be the combined set of the two lists?
                uniprot_ids = locus_subdict["uniprotIds"]
            else:
                uniprot_ids = rev_uniprot_ids
        
            # loop over uniprot IDs found via reverse lookup or during parsing
            for id_ in uniprot_ids:
                # append to output_file
                with open(output_file, "a") as out_tab:
                    out_tab.write(f"{self.ENA_ID}\t{id_}\t{locus}\t{self.CHR}\t{locus_subdict['DIR']}\t{locus_subdict['START']}\t{locus_subdict['END']}\n")


def process_file(
        file_path: str, 
        database_connection,
        output_file: str):
    """
    read gzipped GenBank or EMBL file, loop over entries, search for lines 
    specifying chromosomes via the ID lines, then gather information for each
    gene locus. Check if a gene lcous' proteinId(s) have been already mapped to
    a UniprotID via a database reverse lookup. Write to output_file the
    organization of the chromosome's gene loci.

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
        r"(?:^FT   CDS\s+)"    # search for "FT   CDS" lines
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

    # create the first Record object that will be empty
    enaRecord = Record(ID = "", CHR = 0)
    
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
            
            # regex search found a new ENA ID line
            if search_results[0] and search_results[1]:
                
                # need to handle the previous Record's data
                enaRecord.process_record(
                    count,
                    database_connection,
                    output_file
                )
                    
                # now handle the new line's data
                # search_results[0] is the ENA Accession ID
                ID = search_results[0]
                
                # search_results[1] is the type of chromosome structure, lienar 
                # or circular; there are some non-standard structures
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
                
                # create the new Record object, currently empty except for the 
                # ID and CHR fields
                enaRecord = Record(ID = ID, CHR = CHR)
                    
            # regex search found a protein_id line
            elif search_results[2]:
                enaRecord.proteinIds.add(search_results[2])
                enaRecord.current_locus['proteinIds'].add(search_results[2])

            # regex search found a UniProtKB line
            elif search_results[3]:
                enaRecord.uniprotIds.add(search_results[3])
                enaRecord.current_locus['uniprotIds'].add(search_results[3])
            
            # regex search matched but tuple is filled with empties, must 
            # have found a "FT   CDS" line
            else:
                # if the current_locus' START and END values are non-zero and 
                # the proteinIds set is not empty, then the current_locus dict 
                # is filled with a previous CDS line's data. Stash this locus' 
                # results into the Record's loci_dict
                if (enaRecord.current_locus["START"] 
                        and enaRecord.current_locus["END"]
                        and enaRecord.current_locus["proteinIds"]):
                    # add the current_locus dict to the full loci_dict, using
                    # the count integer as the key
                    enaRecord.add_locus(count)
                    # add one to the count since we're moving on to the next 
                    # locus
                    count += 1
                
                # for the new locus:
                # parse the CDS line to see if it fits the expected regex format
                # find the start and stop values for the sequence as well as the
                # directionality
                cds_result = cds_pattern.findall(line)
                if cds_result:
                    # the two regex groups
                    START, END = cds_result[0]
                    # determine the directionality of the CDS
                    DIR = 0 if "complement" in line else 1
                    # add info to the Record's current_locus dict
                    enaRecord.current_locus["DIR"] = DIR
                    enaRecord.current_locus["START"] = int(START)
                    enaRecord.current_locus["END"] = int(END)
                
    # done parsing file, but haven't finished parsing the last Record object
    enaRecord.process_record(
        count,
        database_connection,
        output_file
    )
    
    return output_file


