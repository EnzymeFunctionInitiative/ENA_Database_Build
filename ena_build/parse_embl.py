
import re
import gzip


###############################################################################
# Global variables
###############################################################################

# create the "ID " search pattern.
# search for ID lines, where group 0 maps to the ENA chromosome ID, group 1 
# maps to the type of chromosome structure (circular or linear or nonstandard 
# strings like XXX), and group 2 maps to the chromosome entry's length in 
# nucleotide base pairs .
# This is a very rigid pattern, creating a list of a tuple of len 3.
ENA_ID_PATTERN = re.compile(r"^ID\s+(\w+);\s\w+\s\w+;\s(\w+);.*;\s(\d+)\sBP") 

# much less rigid, but creates a list of two tuples of len 3 each. 
#ENA_ID_PATTERN = re.compile(r"^ID\s+(\w+);\s\w+\s\w+;\s(\w+);|(\d+)\sBP") 
## original ID line regex pattern, creating a list of a tuple of len 2.
#ENA_ID_PATTERN = re.compile(r"^ID\s+(\w+);\s\w+\s\w+;\s(\w+);") 

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

# If the line matches this search string, then a list of tuples with integer 
# substrings will be returned. Each tuple element in the list is one of the 
# location ranges for the associated ID line. Else, the regex search will 
# return an empty list. 
# See https://www.insdc.org/submitting-standards/feature-table/#3.4 for 
# documentation on Location descriptors and operators. 
## NOTE: The below pattern does not parse e.g. `102`, `102.110`, `102^112`, 
# or any combination of these "ranges" in location substrings. These are 
# uninteresting for our purposes but it is important to note. 
CDS_LOC_PATTERN = re.compile(r"(\d+)\.\.\>?(\d+)")

## original regex pattern:
#CDS_LOC_PATTERN = re.compile(r"(\d+)\..*\.\>?(\d+)")
## grabs the outer two numbers separated by any characters
## ChatGPT was used to refine the below regex pattern. 
#CDS_LOC_PATTERN = re.compile(r"(?<![A-Za-z0-9]\.)\b\d+\b")
## Uses a negative lookbehind for alphanumeric character followed by a period.
## The \b\d+\b captures any whole "word" numbers, reports these matches in a 
## list of ints

# If the line matches this search string, this indicates the start of a new 
# feature block. There are a large number of lines that could be identified
FT_START_PATTERN = re.compile(r"^FT\s\s\s[a-zA-Z0-9-]")

###############################################################################
# Define the Record Object
###############################################################################

class Record():
    def __init__(
            self, 
            ena_id: str, 
            chr_struct: int, 
            chr_len: int, 
            file_path: str) -> None:
        """ 
        Instantiate a Record object for a new chromosome block in an EMBL flat
        file.

        Parameters
        ----------
            ena_id
                str, EMBL/ENA accession ID for the chromosome.
            chr_struct
                int, 0 if chromosome structure/type is "linear" or 1 if 
                "circular".
            chr_len
                int, length of the chromosome as reported on the ID line
            file_path
                str, file path within which the ENA ID entry originates

        Attributes
        ----------
            file_path
                str, file path within which the ENA ID entry originates
            ena_id
                str, EMBL/ENA accession ID for the chromosome.
            chr_struct
                int, 0 if chromosome structure/type is "linear" or 1 if 
                "circular".
            count
                incremented count attribute to denote CDS locus position in the
                chromosome.
            uniprotIds
                set of unique UniProt Accession Ids associated with CDS loci.
            proteinIds
                set of unique protein Ids associated with CDS loci.
            loci_dict
                dict that will contain information for all CDS loci associated
                with the chromosome.
                - keys in loci_dict will be the locus' count value with values
                  being the associated Locus object with "direction", "start",
                  "end", "uniprotIds", and "proteinIds" instance attributes.
            current_locus_lines
                list, used to gather all lines associated with the 
                currently-being parsed CDS block.

        """
        self.file_path = file_path
        self.ena_id = ena_id
        self.chr_struct = chr_struct
        self.length = chr_len
        self.count = 1
        self.uniprotIds = set()
        self.proteinIds = set()
        self.loci_dict = {}
        self.current_locus_lines = []

    def add_locus(self) -> None:
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
        cds_line = "".join(self.current_locus_lines).split("/")[0]
        # remove "FT ", "CDS " and white space from the cds_line string
        for substring in ["FT ","CDS ","\n"," "]:
            cds_line = cds_line.replace(substring,"")
        # feed the cds_line into the CDS_LOC_PATTERN to get a list of tuples of
        # strs associated with the "x..y" formatted ranges. 
        loc_ranges = CDS_LOC_PATTERN.findall(cds_line)
        # make sure the pattern was matched
        if loc_ranges:
            # convert the cds_result list of tuple of strs into list of tuple
            # of ints
            loc_ranges = [tuple(int(i) for i in tup) for tup in loc_ranges]
            # assign the matched groups to the Locus object's attributes
            start, end = process_location_ranges(
                loc_ranges, 
                self.chr_struct, 
                self.length
            )
            direction = 0 if "complement" in cds_line else 1
        # if the cds_line fails to parse, end the method early and clean the
        # locus lines list.
        else:
            print("!!! FT CDS line block failed to be processed. "
                    + f"{self.file_path}:\n{self.current_locus_lines}")
            self.current_locus_lines = []
            return

        # prep some containers for the parsed results
        uniprotIds = set()
        proteinIds = set()

        # Loop over each line in the input locus_lines to find the important
        # optional qualifier lines for UniprotId and proteinId.
        for line in self.current_locus_lines:
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

        # add the Locus object to the Record's loci_dict, using the
        # Record.count attribute as the key.
        if self.count not in self.loci_dict.keys():
            self.loci_dict[self.count] = Locus(
                direction = direction,
                start = start,
                end = end,
                uniprotIds = uniprotIds,
                proteinIds = proteinIds
            )
            # increment the counter for the next Locus to be found
            self.count += 1
        
        # refresh the current_locus_lines string.
        self.current_locus_lines = []
        return

    def process_record(self, db_cnx, output_file) -> None:
        """
        Given a Record object, do one final check for adding a locus then 
        process the loci in the Record object, writing info associated with 
        any loci that are associated with a UniProtId out to a tab separated
        file.

        Parameters
        ----------
            db_cnx
                mysql_database.IDMapper connection object
            output_file
                string or pathlib.Path, file to be written with results from the
                function
        
        """
        # before doing any processing, check to make sure self.ena_id is not an
        # empty string; empty string for self.ena_id is used to denote Records
        # that should not be processed.
        if not self.check_ena_id():
            return

        # check to see if the Record.current_locus_lines is not empty
        if self.check_locus_lines():
            self.add_locus()
        
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
        for locus_count, locus in self.loci_dict.items():
            # gather uniprotIds from the reverse_mapping call for each 
            # proteinId associated with the locus; 
            rev_uniprot_ids = [id_ for proteinId in locus.proteinIds for id_ in ids_mapping.get(proteinId,[]) if proteinId not in no_match]
            ## equivalent to:
            #rev_uniprot_ids = []
            #for proteinId in locus_subdict["proteinIds"]:
            #   for id_ in ids_mapping.get(proteinId,[]):
            #       if proteinId not in no_match:
            #            rev_uniprot_ids.append(id_)

            # check whether the rev_uniprot_ids list is empty
            if not rev_uniprot_ids:
                # if it is empty, use the loci's uniprotIds value instead
                uniprot_ids = locus.uniprotIds
            else:
                uniprot_ids = rev_uniprot_ids
        
            # loop over uniprot IDs found via reverse lookup or during parsing
            for id_ in uniprot_ids:
                # append to output_file
                with open(output_file, "a") as out_tab:
                    out_tab.write(f"{self.ena_id}\t{id_}\t{locus_count}\t{self.chr_struct}\t{locus.direction}\t{locus.start}\t{locus.end}\n")

    def check_ena_id(self):
        """ Check whether the ena_id attribute is an empty string or not """
        return bool(self.ena_id)

    def check_locus_lines(self):
        """ Check whether the current_locus_lines is an empty list or not """
        return bool(self.current_locus_lines)


###############################################################################
# Define the Locus Object
###############################################################################

class Locus():
    def __init__(self, direction: int, start: int, end: int, uniprotIds: set, 
                 proteinIds: set) -> None:
        """
        Instantiate a Locus object used to gather associated Locus data.

        Parameters
        ----------
            direction
                int, integer representation of the direction of the CDS on the
                antiparallel strands of DNA. Direction is relative to the 
                presented strand for the full ENA ID Record. Values of 0 if on
                the complement stand, else 1.
            start
                int, position of starting nucleic acid base associated with the
                CDS Locus.
            end
                int, position of the ending nucleic acid base associated with
                the CDS Locus.
            uniprotIds
                set of strings, UniProtKB accession IDs. 
            proteinIds
                set of strings, protein identifier IDs, including the version 
                number. 
        
        Attributes
        ----------
        See above.
        """
        self.direction = direction
        self.start = start
        self.end = end
        self.uniprotIds = uniprotIds
        self.proteinIds = proteinIds

###############################################################################
# Parsing an EMBL flat file
###############################################################################

def process_id_line(line: str, file_path: str) -> tuple[str, int]:
    """
    Given the line, parse it as if it is an "ID" line from the embl flat file.

    Parameters
    ----------
        line
            str, the full line string that starts with "ID".
        file_path
            str, file path within which the "ID" line originates.

    Returns
    -------
        ena_id
            str, the EMBL Accession ID for the chromosome entry being parsed.
        chr_struct
            int, 0 if chromosome structure/type is "linear" or 1 if "circular".
        chr_len
            int, length in nucleotide base pairs of the chromosome entry.
    """
    # parse ID line using regex to get list of tuple of len 3, group 0 is the 
    # ENA ID, group 1 is the chromosome structure type, group 2 is the length
    # if the pattern is not matched, id_groups is an empty list.
    id_groups = ENA_ID_PATTERN.findall(line)
    if id_groups:
        ena_id, chr_struct, chr_len = id_groups[0]
        
        # chr_struct is the type of chromosome structure, linear or circular; there 
        # are some non-standard structures, skip those.
        if chr_struct in ["linear","circular"]:
            # check if the type is either linear or circular
            chr_struct = 1 if chr_struct == "linear" else 0
        else:
            print("!!! Unknown chromosome type observed in "
                    + f"{file_path}:{line.strip()}")
            # by replacing ena_id with an empty string, we are effectively 
            # ignoring unexpected chromosome type strings; we'll ignore the 
            # "" ena_id string in the Record.process_record() method
            ena_id = ""
            chr_struct = -1
            chr_len = 0

    else:
        # If the ENA_ID_PATTERN isn't matched then we've hit an unexpected 
        # result. Log it and return uninteresting values. 
        ena_id = ""
        chr_struct = -1
        chr_len = 0
        print("!!! Ill-formatted ID line observed in " 
                + f"{file_path}:{line.strip()}")

    return ena_id, chr_struct, int(chr_len)


def process_location_ranges(
        loc_ranges: list[tuple[int,int]],
        linear_chromosome: bool,
        chromosome_length: int) -> list[int,int]:
    """
    Given a list of sequence position ranges, determine the real start and end
    positions for the associated locus. This takes into account the chromosome
    structure and length since circular chromosome structures complicate the 
    matter.

    Parameters
    ----------
        loc_ranges
            list of tuples of ints, each tuple is a pair of positions. 
        linear_chromosome
            boolean or equivalent, signifies if the ranges were pulled from a 
            linear chromosome structure or not. 
        chromosome_length
            int, length of the full chromosome; used to determine if a range 
            spans the relative start position of a circular chromosome

    Returns
    -------
        start, end
            ints, the appropriate start and stop indices for the given range(s)
            and chromosome structure type and length.
    """
    if linear_chromosome:
        # handle the ranges assuming the hard boundary conditions of 1 and 
        # chromosome_length
        flattened_ranges = sum(loc_ranges,())
        return min(flattened_ranges), max(flattened_ranges)
    else:
        # note: sequence position indices are 1-indexed rather than 0-indexed.
        # loop over all ranges to see if the chromosome length is in 
        # one of the ranges.
        passes_length = [
            chromosome_length in range(elem[0],elem[1]+1) 
            for elem in loc_ranges
        ]
        # loop over all ranges to see if the start index (1) is in one of the 
        # ranges.
        includes_start = [
            1 in range(elem[0],elem[1]+1) 
            for elem in loc_ranges
        ]
        if True in passes_length and True in includes_start:
            # found an instance where the locus spans the chromosome length. 
            # !!! assume that the join components (required for CDSs that span
            # a relative start position) are in ascending order, crossing the
            # cyclic boundary condition. 
            # NOTE: potential source of bugs... And this is really the important
            # logic to consider for this function.
            # is it possible to have "join(1..27,len-100..len)"? If so,
            # the return value will be wrong since its assuming ascending order
            # of the ranges.
            start = loc_ranges[0][0]
            end = loc_ranges[-1][-1]
            if start < end: 
                print(f"!!! {loc_ranges}, {linear_chromosome}, {chromosome_length}")
            return start, end
        else:
            # no range passes through the chromosome_length so treat like a 
            # linear sequence
            flattened_ranges = sum(loc_ranges,())
            return min(flattened_ranges), max(flattened_ranges)


def process_file(
        file_path: str, 
        database_connection,
        output_file: str) -> str:
    """
    Read gzipped GenBank or EMBL file, loop over entries, search for lines 
    specifying chromosomes via the ID lines, then gather information for each
    gene locus. Check if a gene locus' proteinId(s) have been already mapped to
    a UniprotID via a database reverse_mapping() call. Write to output_file the
    organization of the chromosome's gene loci.

    Parameters
    ----------
        file_path
            str or pathlib.Path, assumed to be a gzipped embl flatfile.
        database_connection
            IDMapper object, connection object to the MYSQL database that 
            will be queried. 
        output_file
            str or pathlib.Path, file to be written with results from the
            function.

    Return
    ------
        output_file
            str or pathlib.Path, file written to with results; if no results
            were found, this file will not actually be written so a check for 
            its existence outside of this function is necessary.
    """
    # create the first Record object that will be empty
    enaRecord = Record(
        ena_id = "", 
        chr_struct = -1, 
        chr_len = 0, 
        file_path = file_path)
    
    # open and read the gzipped file
    with gzip.open(file_path, 'rt') as f:
        # loop over each line in f without reading the whole file
        for line in f: 
            # check for lines that do not start with "FT   ", "ID   ", nor 
            # "OC   ", we currently don't do anything with those lines so just
            # move on
            if not line.startswith(("FT   ", "ID   ", "OC   ")):
                continue

            # check for "ID   " lines to grab ENA IDs and chromosome structure
            # type and len. Before doing that, need to check whether a previous
            # Record (and associated Locus) is ready to be processed. 
            elif line.startswith("ID   "):
                # check if the enaRecord.current_locus_lines is full
                if enaRecord.check_locus_lines():
                    # parse the string and add the locus to the loci_dict
                    enaRecord.add_locus()

                # check if the enaRecord object was filled with info before
                # processing it
                if enaRecord.check_ena_id():
                    # process the previous Record's data
                    enaRecord.process_record(
                        database_connection,
                        output_file
                    )
                
                # process the current ID line to grab the ena_id and 
                # chr_struct info
                ena_id, chr_struct, chr_len = process_id_line(line, file_path)

                # create the new Record object that will start out empty except
                # for the ena_id, chr_struct, and chr_len instance attributes
                enaRecord = Record(
                    ena_id = ena_id,
                    chr_struct = chr_struct,
                    chr_len = chr_len,
                    file_path = file_path
                )

            # only interested in parsing Records associated with the Fungi 
            # kingdom when considering genomes from Eukaryota domain. 
            # an OC line is always found after an ID line and before any FT 
            # lines. Just overwrite the enaRecord object if a non-Fungi
            # eukaryote is being parsed.
            elif (line.startswith("OC   ") 
                    and "Eukaryota" in line 
                    and " Fungi" not in line):
                enaRecord = Record(
                    ena_id = "",
                    chr_struct = -1,
                    chr_len = 0,
                    file_path = ""
                )

            # check whether the current enaRecord object has a True-like ena_id
            # attribute. If it doesn't, then any lines can be skipped since the
            # active Record will not be parsed.
            elif not enaRecord.check_ena_id():
                continue

            # check for the start of _any_ chromosome feature block. This 
            # includes a "FT   CDS " line.
            elif FT_START_PATTERN.match(line):
                # since we've hit the start of a new feature block, we need to
                # check that the previous block was an "FT   CDS" block. Do so
                # by checking whether the Record.current_locus_lines is an
                # empty list. If it isn't, then the lines for a new locus have
                # been fully gathered and need to be processed into a Locus
                # object.
                if enaRecord.check_locus_lines():
                    # process the CDS block string
                    enaRecord.add_locus()
                
                # check for new gene coding sequence blocks
                if line.startswith("FT   CDS "):
                    # start the Record.current_locus_lines string
                    enaRecord.current_locus_lines.append(line)
            
            # At this point, the only lines remaining are "FT\s+" lines that 
            # may be associated with a CDS block or not. Need to check if the 
            # Record.current_locus_lines is not empty.
            elif enaRecord.check_locus_lines() and line.startswith("FT    "):
                enaRecord.current_locus_lines.append(line)
 
    # done parsing file, but haven't finished parsing the last Record object
    enaRecord.process_record(database_connection, output_file)
    
    return output_file


