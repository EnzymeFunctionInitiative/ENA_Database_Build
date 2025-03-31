
# Tests:
- `regex_test.py`, runs small scale testing of the four regex patterns used in `../ena_build/parse_embl.py`. 
- `location_parsing_test.py`, runs the location regex and parsing functions for a small but complete set of CDS location strings. Does so for both linear and circular chromosome structure types since the logic for parsing the location substrings is different. 

# To be developed:
- Unit tests for creating, adding to, and processing Record objects. * 
- Unit tests for creating and handling Locus objects. 
- Unit tests for creating, executing queries, and closing an IDMapper object. *
- Unit tests for the glob style task functions defined in `dask_tasks.py`. 
- Integration tests for parsing a real EMBL gzipped flat file (dat.gz) via the `process_file()` function. * 
- Integration tests for the full `dask_tskmgr.py` workflow. *

"*" denotes tests that will require a connection to the EFI-DB MySQL server, so won't be easily implementable for outside parties. 

