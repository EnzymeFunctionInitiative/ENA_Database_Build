
import re
import gzip


###############################################################################
# Global variables
###############################################################################

# create the "ID " search pattern.
# search for ID lines, group 0 and 1 map to ID and type of chromosome
# structure (circular or linear or nonstandard strings like XXX);
ENA_ID_PATTERN = re.compile(r"^ID\s+(\w+);\s\w+\s\w+;\s(\w+);") 

# search for "protein_id" FT lines, group 0 maps to the foreign ID string. 
# search for UniProtKB/... database accession ID lines, group 1 matches the
# associated accession ID. 
XREF_SEARCH_STRS = [
    r'^FT\s+\/protein_id=\"([a-zA-Z0-9\.]+)\"', 
    r'^FT\s+\/db_xref=\"UniProtKB\/[a-zA-Z0-9-]+:(\w+)\"'
]
# compile the combined search pattern.
# for each line that successfully matches one of the search strings, a list
# of one tuple of len 2 will be created. The zeroth element of the tuple
# maps to a protein_id string and the first maps to a UniProtKB accession 
# ID. Depending on which line is matched, one of these elements 
# will be an empty string _but_ the tuple will always be len 2. If the line 
# does not match either search strings, then an empty list will be returned. 
XREF_SEARCH_PATTERN = re.compile("|".join(XREF_SEARCH_STRS))

# If the line matches this search string, then the a list of a tuple of 
# len 2 will be returned. The zeroth and first elements will be the START
# and END values for the sequence, respectively. Else, the regex search 
# will return an empty list. 
CDS_LOC_PATTERN = re.compile(r"(\d+)\..*\.\>?(\d+)")

# If the line matches this search string, this indicates the start of a new 
# feature block. There are a large number of lines that could be identified
FT_START_PATTERN = re.compile(r"^FT\s\s\s[a-zA-Z0-9-]")

###############################################################################
# Define the Record Object
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
            loci_dict
                dict that will contain information for all CDS loci associated
                with the chromosome.
                - keys in loci_dict will be the locus' count value with values
                  being a subdict with "DIR", "START", "END", "uniprotIds", and
                  "proteinIds" keys and associated values for the locus.
            current_locus_lines
                str, used to gather all lines associated with the 
                currently-being parsed CDS block.

        """
        self.ENA_ID = ID
        self.CHR = CHR
        self.count = 1
        self.uniprotIds = set()
        self.proteinIds = set()
        self.loci_dict = {}
        self.current_locus_lines = ""

    def add_locus(self):
        """
        Parse the self.current_locus_lines string and add the processed 
        results to the Record's loci_dict subdict. Refresh the 
        self.current_locus_lines string at the end.
        """
        # We need to parse the first or more lines that contain the location
        # information for the coding sequence. The main separater between the
        # end of the "FT   CDS" line(s) and the next line (whatever that line
        # may be) is a forward slash, which denotes optional qualifiers. See
        # https://www.ebi.ac.uk/ena/WebFeat/ and 
        # https://www.insdc.org/submitting-standards/feature-table/#3.3 for
        # more details.
        # this feels pretty quick and dirty to me...
        cds_line = self.current_locus_lines.split('/')[0]
        # there's a whole bunch of annoying "FT ", "CDS " and white space
        # elements in the cds_line, so get rid of those.
        for substring in ["FT ","CDS ","\n"," "]:
            cds_line = cds_line.replace(substring,"")
        # feed the cleaned up version of the cds_line into the CDS_LOC_PATTERN
        cds_result = CDS_LOC_PATTERN.findall(cds_line)
        # make sure the pattern was matched
        if cds_result:
            # assign the matched groups to the Locus object's attributes
            START, END = cds_result[0]
            DIR = 0 if "complement" in cds_line else 1
        # if the cds_line fails to parse, end the method early.
        else:
            self.current_locus_lines = ""
            print(f"!!! FT CDS line block failed to be processed. {file_path}:"
                    + f"\n{self.current_locus_lines}")
            return 1
            #temp = self.current_locus_lines
            #raise InvalidLocusException(temp)

        # prep some containers for the parsed results
        uniprotIds = set()
        proteinIds = set()

        # Loop over each line in the input locus_lines to find the important
        # optional qualifier lines for UniprotId and proteinId.
        for line in self.current_locus_lines.split("\n"):
            # throw the line into the XREF_SEARCH_PATTERN
            search_results = XREF_SEARCH_PATTERN.findall(line)
            # make sure the pattern was matched, otherwise move on
            if search_results:
                # if the regex pattern is matched, search_results will be a
                # list of a tuple of len 2. One or the other element in this
                # tuple will be an empty string since no FT lines contain both
                # uniprot or protein IDs.
                proteinId, uniprotId = search_results[0]
                if uniprotId:
                    uniprotIds.add(uniprotId)
                    self.uniprotIds.add(uniprotId)
                elif proteinId:
                    proteinIds.add(proteinId)
                    self.proteinIds.add(proteinId)

        # add the Locus to the Record's loci_dict, using the Record.count
        # attribute as the key.
        if self.count not in self.loci_dict.keys():
            self.loci_dict[self.count] = {
                "uniprotIds": uniprotIds,
                "proteinIds": proteinIds,
                "DIR": DIR,
                "START": START,
                "END": END
            }
            # increment the counter for the next Locus to be found
            self.count += 1
        
        # refresh the current_locus_lines string.
        self.current_locus_lines = ""
        return 0

    def process_record(self, db_cnx, output_file):
        """
        given a Record object, do one final check for adding a locus then 
        process the loci in the Record object, writing info associated with 
        any loci that are associated with a UniProtId out to a tab separated
        file

        PARAMETERS
        ----------
            db_cnx
                mysql_database.IDMapper connection object
            output_file
                string or pathlib.Path, file to be written with results from the
                function
        
        """
        # before doing any processing, check to make sure self.ENA_ID is not an
        # empty string; empty string for self.ENA_ID is used to denote Records
        # that should not be processed.
        if not self.ENA_ID:
            return

        # check to see if the Record.count is not a key in the Record's 
        # loci_dict already 
        if (self.count not in self.loci_dict.keys()
            rt = self.add_locus()
            ## check for non-zero return code
            #if rt != 0:
            #    # do something
        
        # use self.proteinIds set as input to the MySQL Database query
        # ids_mapping is a dict of proteinIds as keys with uniprotIds as values
        # no_match is a list of proteinIds that did not map to a uniprotId in 
        # the database. 
        ids_mapping, no_match = db_cnx.reverse_mapping(
            list(self.proteinIds)
        )

        # loop over each locus in the Record object, get the mapping btw 
        # the proteinIds and uniprotIds from the reverse_mapping results
        # Write to file if the locus has an associated uniprotId.
        for locus, locus_subdict in self.loci_dict.items():
            # gather uniprotIds from the reverse_mapping call for each 
            # proteinId associated with the locus; 
            rev_uniprot_ids = [id_ for proteinId in locus_subdict["proteinIds"] for id_ in ids_mapping.get(proteinId,[]) if proteinId not in no_match]
            ## equivalent to:
            #rev_uniprot_ids = []
            #for proteinId in locus_subdict["proteinIds"]:
            #    if proteinId not in no_match:
            #        for id_ in ids_mapping.get(proteinId,[]):
            #            rev_uniprot_ids.append(id_)

            # check whether the rev_uniprot_ids list is empty
            if not rev_uniprot_ids:
                # if it is empty, use the loci's uniprotIds value instead
                uniprot_ids = locus_subdict["uniprotIds"]
            else:
                uniprot_ids = rev_uniprot_ids
        
            # loop over uniprot IDs found via reverse lookup or during parsing
            for id_ in uniprot_ids:
                # append to output_file
                with open(output_file, "a") as out_tab:
                    out_tab.write(f"{self.ENA_ID}\t{id_}\t{locus}\t{self.CHR}\t{locus_subdict['DIR']}\t{locus_subdict['START']}\t{locus_subdict['END']}\n")

            # remove the locus key from the Record.loci_dict to keep this clean
            self.loci_dict.pop(locus)


###############################################################################
# Parsing an EMBL flat file
###############################################################################

def process_id_line(line: str):
    """
    """
    # parse ID line using regex to get list of tuple of len 2 group 0 is the 
    # ENA ID, group 1 is the chromosome structure type
    id_groups = ENA_ID_PATTERN.findall(line)
    if id_groups:
        ID, CHR = id_groups[0]
        
        # CHR is the type of chromosome structure, linear or circular; there 
        # are some non-standard structures, skip those.
        if CHR in ["linear","circular"]:
            # check if the type is either linear or circular
            CHR = 1 if CHR == "linear" else 0
        else:
            print("!!! Unknown chromosome type observed in "
                    + f"{file_path}:{line}")
            # by replacing ID with an empty string, we are effectively 
            # ignoring unexpected chromosome type strings; we'll ignore the 
            # "" ID string in the Record.process_record() method
            ID = ""
            CHR = -1
    else:
        # NOTE: NEED TO CHECK WHICH LINES ARE BEING CAUGHT HERE
        # DO THEY MAP TO AN ACTUAL ID LINE THAT DOESN'T MAP TO THE 
        # ID_PATTERN? DO THEY SIGNIFY THE START OF A NEW RECORD 
        # OBJECT OR CAN I JUST `CONTINUE` TO THE NEXT RECORD
        ID = ""
        CHR = -1
        print("!!! Ill-formatted ID line observed in " 
                + f"{file_path}:{line}")

    return ID, CHR


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
    # create the first Record object that will be empty
    enaRecord = Record(ID = "", CHR = 0)
    
    # open and read the gzipped file
    with gzip.open(file_path, 'rt') as f:
        # loop over each line in f without reading the whole file
        for line in f: 
            # check for lines that do not start with "FT   ", "ID   ", nor 
            # "OC   ", we currently don't do anything with those lines so just
            # move on
            if not line.startswith(("FT   ", "ID   ", "OC   ")):
                continue

            # check for "ID   " lines to grab ENA IDs and Chromosome structure
            # type
            elif line.startswith("ID   "):
                # check if the enaRecord.current_locus_lines is full
                if enaRecord.current_locus_lines:
                    # parse the string and add the locus to the loci_dict
                    rt = enaRecord.add_locus()
                    ## check for non-zero return code
                    #if rt != 0:
                    #    # do something

                # check if the enaRecord object was filled with info before
                # processing it
                if not enaRecord.ENA_ID:
                    # need to handle the previous Record's data
                    enaRecord.process_record(
                        database_connection,
                        output_file
                    )
                
                # process the ENA ID line to grab the ID and CHR info
                ID, CHR = process_id_line(line)

                # create the new Record object, currently empty except for the 
                # ID and CHR fields
                enaRecord = Record(ID = ID, CHR = CHR)
                continue

            # only interested in parsing Records associated with the Fungi 
            # kingdom when considering genomes from Eukaryota domain. 
            # an OC line is always found after an ID line and before any FT 
            # lines. Just overwrite the enaRecord object.
            elif (line.startswith("OC   ") 
                    and "Eukaryota" in line 
                    and " Fungi" not in line):
                #print(f"!!! Found non-fungi eukaryote in {file_path}: {enaRecord.ENA_ID}, {line}")
                enaRecord = Record(ID = "", CHR = -1)
                continue

            # check whether the current enaRecord object has a True-like ENA_ID
            # attribute. If it doesn't, then any lines can be skipped since the
            # active Record will not be parsed into a file.
            elif not enaRecord.ENA_ID:
                continue

            # check for the start of _any_ chromosome feature block. This 
            # includes a "FT   CDS " line.
            elif FT_START_PATTERN.match(line):
                # since we've hit the start of a new feature block, we need to
                # check that the previous block was an "FT   CDS" block. Do so
                # by checking whether the Record.current_locus_lines is an
                # empty string. If it isn't, then the lines for a new locus
                # have been fully gathered and need to be processed into a 
                # Locus object.
                if enaRecord.current_locus_lines:
                    # process the CDS block string
                    rt = enaRecord.add_locus()
                    ## check for non-zero return code
                    #if rt != 0:
                    #    # do something
                
                # check for new gene coding sequence blocks
                if line.startswith("FT   CDS "):
                    # start the Record.current_locus_lines string
                    enaRecord.current_locus_lines = line
                    continue
            
            # At this point, the only lines remaining are "FT\s+" lines that 
            # may be associated with a CDS block or not. Need to check if the 
            # Record.current_locus_lines is not empty.
            elif enaRecord.current_locus_lines and line.startswith("FT    "):
                enaRecord.current_locus_lines += line
 
    # done parsing file, but haven't finished parsing the last Record object
    enaRecord.process_record(database_connection, output_file)
    
    return output_file


