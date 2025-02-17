
import re
import gzip

###############################################################################
# Parsing an EMBL flat file
###############################################################################

class Record():
    def __init__(self, ID: str, CHR: int):
        """ 
        Instantiate a Record object for a new chromosome block in an EMBL flat 
        file

        PARAMETERS
        ----------
            ID
                str, EMBL/ENA accession ID for the chromosome
            CHR
                int, 0 if chromosome structure/type is "linear" or 1 if 
                "circular"

        ATTRIBUTES
        ----------
            ID
                str, EMBL/ENA accession ID for the chromosome
            CHR
                int, 0 if chromosome structure/type is "linear" or 1 if 
                "circular"
            count
                incremented count attribute to denote CDS locus position in the
                chromosome
            uniprotIds
                set of unique UniProt Accession Ids associated with CDS loci 
            proteinIds
                set of unique "foreign" protein Ids associated with CDS loci
            loci_dict and current_locus
                dicts that will contain information for all and current CDS
                loci associated with the chromosome, respectively.
                - keys in loci_dict will be the locus' count value with values
                  being a subdict with "DIR", "START", "END", "uniprotIds", and
                  "proteinIds" keys and associated values for the locus.

        """
        self.ENA_ID = ID
        self.CHR: CHR,
        self.count: 0,
        self.uniprotIds: set(),
        self.proteinIds: set(),
        # these two dicts will be filled with CDS entries, keys being the count
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

    def add_locus(self, locus_id):
        """
        Add the current locus' dict to the chromosome's loci_dict then refresh
        the current_locus dict
        
        PARAMETERS
        ----------
            locus_id
                str or int, used as the key in self.loci_dict that maps to the
                current_locus subdict.
        """
        self.loci_dict[locus_id] = self.current_locus
        self.current_locus = {
            "uniprotIds": set(),
            "proteinIds": set(),
            "DIR": -9999,
            "START": 0,
            "END": 0,
        }

    def process_record(self, last_locus, db_cnx, output_file):
        """
        given a Record object, do one final check for adding a locus then 
        process the loci in the Record object, writing info associated with 
        any loci that are associated with a UniProtId out to a tab separated
        file

        PARAMETERS
        ----------
            last_locus
                str or int, used as the key in self.loci_dict that maps to the
                current_locus subdict.
            db_cnx
                mysql_database.IDMapper connection object
            output_file
                string or pathlib.Path, file to be written with results from the
                function
        
        """
        # check to see if the last_locus is not a key in the Record's loci_dict
        # as well as that locus' START, END, and DIR values are different from
        # the default/starting values
        if (last_locus not in self.loci_dict.keys()
                and self.current_locus["START"]
                and self.current_locus["END"]
                and self.current_locus["DIR"] >= 0):
            self.add_locus(last_locus)
        
        # loop over each locus in the Record object, doing a reverse_lookup on
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
    # create the regex search strings before hand to improve efficiency of
    # doing the regex searches on lines of the file(s).
    # search for "protein_id" FT lines, group 0 maps to the foreign ID string. 
    # search for UniProtKB/... database accession ID lines, group 1 matches the
    # associated accession ID. 
    search_strs = [
        r'(?:^FT\s+\/protein_id=\"([a-zA-Z0-9\.]+)\")', 
        r'(?:^FT\s+\/db_xref=\"UniProtKB\/[a-zA-Z0-9-]+:(\w+)\")'
    ]
    
    # compile the combined search pattern.
    # for each line that successfully matches one of the search strings, a list
    # of one tuple of len 2 will be created. The zeroth element of the tuple
    # maps to a protein_id string and the first maps to a UniProtKB accession 
    # ID. Depending on which line is matched, some or all of these elements 
    # could be empty strings _but_ the tuple will always be len 2. If the line 
    # does not match any search strings, then an empty list will be returned. 
    search_pattern = re.compile("|".join(search_strs))
  
    # create the "ID " search pattern.
    # search for ID lines, group 0 and 1 map to ID and type of chromosome
    # structure (circular or linear or nonstandard strings like XXX);
    id_pattern = re.compile(r"(?:^ID\s+(\w+);\s\w\w\s\w;\s(\w+);\s.*)") 

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
            # check for lines that do not start with "FT   " nor "ID   "
            # we currently don't do anything with those lines so just move on
            if not line.startswith(("FT   ", "ID   ")):
                continue

            # check for "ID" lines
            elif line.startswith("ID   "):
                
                # need to handle the previous Record's data
                enaRecord.process_record(
                    count,
                    database_connection,
                    output_file
                )
                
                # parse ID line using regex to get groups
                # group 0 is the ENA ID, group 1 is the chromosome structure 
                # type
                ID, CHR = id_pattern.findall(line)[0]

                # CHR is the type of chromosome structure, linear or circular; 
                # there are some non-standard structures, skip those.
                if CHR in ["linear","circular"]:
                    # check if the type is either linear or circular
                    CHR = 1 if CHR == "linear" else 0
                else:
                    print(f"!!! Unknown chromosome type observed in {file_path}: {line}")
                    # by replacing ID with an empty string, we are effectively 
                    # ignoring unexpected chromosome type strings; we'll
                    # ignore the "" ID string in the Record.process_record()
                    # method
                    ID = ""
                
                # create the new Record object, currently empty except for the 
                # ID and CHR fields
                enaRecord = Record(ID = ID, CHR = CHR)
                continue

            # check for new gene coding sequence blocks
            elif line.startswith("FT   CDS "):
                # if the count value is not already a key in the 
                # Record.loci_dict and the current_locus' START and END values 
                # are non-zero and the DIR value is not the default, then the 
                # current_locus dict is filled with a previous CDS line's data. 
                # Stash this locus' results into the Record's loci_dict.
                if (count not in enaRecord.loci_dict.keys()
                        and enaRecord.current_locus["START"]
                        and enaRecord.current_locus["END"]
                        and enaRecord.current_locus["DIR"] >= 0):
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
                
                continue
            
            # now we need to consider the diverse set of "FT\s+" formatted lines
            # regex still feels like the fastest way to check for the lines that
            # match the search patterns
            search_results = search_pattern.findall(line)
            if search_results:
                uniprotId, proteinId = search_results[0]
                # one or the other will be an empty string since no FT lines
                # contain both uniprot or protein IDs
                if uniprotID:
                    enaRecord.uniprotIds.add(uniprotId)
                    enaRecord.current_locus['uniprotIds'].add(uniprotId)
                elif proteinId:
                    enaRecord.proteinIds.add(proteinId)
                    enaRecord.current_locus['proteinIds'].add(proteinId)
 
    # done parsing file, but haven't finished parsing the last Record object
    enaRecord.process_record(count, database_connection, output_file)
    
    return output_file


