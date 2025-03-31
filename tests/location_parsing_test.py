
import pytest

import parse_embl

def test_linear_loc_parsing():
    """ 
    Testing wrapper func for the location parsing function, given a linear 
    chromosome structure type.
    """
    chr_struct = 1
    chromosome_length = 5163

    # simple example
    line = "FT   CDS             <1..1701"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,1701)

    # completely ignoring the reference entry info
    line = "FT   CDS             join(1..100,J00194.1:100..202)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # singular nucleotide in the join list
    line = "FT   CDS             join(1..100,J00194.1:100..202,300)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # difficult bc join's ranges are unsorted
    line = "FT   CDS             join(complement(4918..5163),complement(2691..4571))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (2691,5163)

#    # m
#    line = ""
#    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
#    matches = [tuple(int(i) for i in tup) for tup in matches]
#    assert parse_embl.process_location_ranges(
#            line, chr_struct, chromosome_length) == (,)


def test_circular_loc_parsing():
    """ 
    Testing wrapper func for the location parsing function, given a circular
    chromosome structure type.
    """
    chr_struct = 0
    chromosome_length = 5163
    
    # simple example
    line = "FT   CDS             <1..1701"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
        matches, chr_struct, chromosome_length) == (1,1701)

    # completely ignoring the reference entry info
    line = "FT   CDS             join(1..100,J00194.1:100..202)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # singular nucleotide in the join list
    line = "FT   CDS             join(1..100,J00194.1:100..202,300)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # difficult bc join's ranges are unsorted, possibly indicating the locus
    # spans the circular chromosome structure's relative start position 
    # BUT! this test should catch on the first range because the 2nd range 
    # doesn't start at 1. 
    # ...
    line = "FT   CDS             join(complement(4918..5163),complement(2691..4571))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (2691,5163)

    # difficult bc join's ranges are unsorted, possibly indicating the locus
    # spans the circular chromosome structure's relative start position 
    # This test should catch the second range because it starts at 1.
    # ...
    line = "FT   CDS             join(complement(4918..5163),complement(1..45))"
    # grab the paired integers from the line
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    # grab the start and stop values from the list of tuples and check against
    # ground truth
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (4918,45)


