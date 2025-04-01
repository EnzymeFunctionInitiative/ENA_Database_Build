
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

    # completely ignoring the reference entry info, may get us into trouble
    line = "FT   CDS             join(1..100,J00194.1:100..202)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # ignore singular nucleotide positions in the join list
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

    # completely ignoring the reference entry info, may get us into trouble
    line = "FT   CDS             join(1..100,J00194.1:100..202)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # ignore singular nucleotide positions in the join list
    line = "FT   CDS             join(1..100,J00194.1:100..202,300)"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,202)

    # difficult bc join's ranges are unsorted, possibly indicating the locus
    # spans the circular chromosome structure's relative start position 
    # BUT! this test does not span the relative start boundary since one range
    # does not include the chr_length, so just take the min and max of the 
    # flattened ranges
    line = "FT   CDS             join(complement(51..60),complement(1..45))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,60)

    # difficult bc join's ranges are unsorted, possibly indicating the locus
    # spans the circular chromosome structure's relative start position 
    # BUT! this test does not span the relative start boundary since one range
    # does not include 1, so just take the min and max of the flattened ranges
    line = "FT   CDS             join(complement(4918..5163),complement(2691..4571))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (2691,5163)

    # difficult bc join's ranges are unsorted, possibly indicating the locus
    # spans the circular chromosome structure's relative start position 
    # Both the chr_len and 1 are found in the list of ranges, indicating the
    # protein does span the chromosome's relative start position.
    line = "FT   CDS             join(complement(4918..5163),complement(1..45))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (4918,45)
    
    ###########################################################################
    # Both the chr_len and 1 are found in the list of ranges, but not in the
    # assumed ascending order that the code currently expects. This is a bug.
    line = "FT   CDS             join(complement(1..45),complement(4918..5163))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,5163) # wrong!

    # Both the chr_len and 1 are found in the list of ranges, but not in the
    # assumed ascending order that the code currently expects. This is a bug.
    line = "FT   CDS             join(complement(1..45),complement(4918..5163),complement(51..60))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (1,60) # wrong!

    # Both the chr_len and 1 are found in the list of ranges, but not in the
    # assumed ascending order that the code currently expects. This is a bug.
    line = "FT   CDS             join(complement(4918..5163),complement(51..60),complement(1..45))"
    matches = parse_embl.CDS_LOC_PATTERN.findall(line)
    matches = [tuple(int(i) for i in tup) for tup in matches]
    assert parse_embl.process_location_ranges(
            matches, chr_struct, chromosome_length) == (4918,45) # wrong!


