
import pytest

import parse_embl

# prep data for all tests
linear_chromosome_struct = 1
circular_chromosome_struct = 0
chromosome_length = 1000

# for the most part, these test examples are fairly self-explanatory for a 
# linear chromosome structure. Just take the max and min of the flattened 
# list to get the start and stop values. The logic is more complex for 
# circular chromosome structures since crossing the periodic boundary can 
# happen. The code detects these cases and correctly parses the potentially
# complex "joins" that may be present.
test_data = [
    # should not see differemt results in the below tests for linear and
    # circular chromosome structures
    pytest.param(
        [(1,700)], linear_chromosome_struct, (1,700), id = "lin_simple"
    ),
    pytest.param(
        [(1,700)], circular_chromosome_struct, (1,700), id = "circ_simple"
    ),
    pytest.param(
        [(1,100),(100,202)], linear_chromosome_struct, (1,202), id = "lin_simple_join"
    ),
    pytest.param(
        [(1,100),(100,202)], circular_chromosome_struct, (1,202), id = "circ_simple_join"
    ),
    pytest.param(
        [(491,516),(269,457)], linear_chromosome_struct, (269,516), id = "lin_unsorted_join"
    ),
    pytest.param(
        [(491,516),(269,457)], circular_chromosome_struct, (269,516), id = "circ_unsorted_join"
    ),
    pytest.param(
        [(491,516),(110,220),(269,457),(518,600)], linear_chromosome_struct, (110,600), id = "lin_many_joins"
    ),
    pytest.param(
        [(491,516),(110,220),(269,457),(518,600)], circular_chromosome_struct, (110,600), id = "circ_many_joins"
    ),
    # edge case, range(s) give full coverage of the chromosome
    pytest.param(
        [(1,1000)], linear_chromosome_struct, (1,1000), id = "lin_full_cov"
    ),
    pytest.param(
        [(1,1000)], circular_chromosome_struct, (1,1000), id = "circ_full_cov"
    ),
    pytest.param(
        [(1,500),(501,1000)], linear_chromosome_struct, (1,1000), id = "lin_join_full_cov"
    ),
    pytest.param(
        [(1,500),(501,1000)], circular_chromosome_struct, (1,1000), id = "circ_join_full_cov"
    ),
    # should see different results in the below tests for linear and circular
    # chromosome structures
    # clearly a gene spanning the boundary
    pytest.param(
        [(1,70),(900,1000)], linear_chromosome_struct, (1,1000), id = "lin_spans_boundary"
    ),
    pytest.param(
        [(1,70),(900,1000)], circular_chromosome_struct, (900,70), id = "circ_spans_boundary"
    ),
    # seq position 1-25 not included in the join ranges... NOT SURE THIS CAN 
    # ACTUALLY HAPPEN FOR GENES SPANNING THE BOUNDARY
    pytest.param(
        [(25,70),(900,1000)], linear_chromosome_struct, (25,1000), id = "lin_skip_start"
    ),
    pytest.param(
        [(25,70),(900,1000)], circular_chromosome_struct, (900,70), id = "circ_skip_start"
    ),
    # seq position 1000 not included in the join ranges... NOT SURE THIS CAN 
    # ACTUALLY HAPPEN FOR GENES SPANNING THE BOUNDARY
    pytest.param(
        [(1,70),(900,999)], linear_chromosome_struct, (1,999), id = "lin_skip_end"
    ),
    pytest.param(
        [(1,70),(900,999)], circular_chromosome_struct, (900,70), id = "circ_skip_end"
    ),
    # both "termini" not included in the join ranges... NOT SURE THIS CAN
    # ACTUALLY HAPPEN FOR GENES SPANNING THE BOUNDARY
    pytest.param(
        [(25,70),(900,999)], linear_chromosome_struct, (25,999), id = "lin_skip_termini"
    ),
    pytest.param(
        [(25,70),(900,999)], circular_chromosome_struct, (900,70), id = "circ_skip_termini"
    ),
    # edge case: near full coverage, suggesting a gene that does not span the 
    # boundary; but since the gap is larger than the wrap_gap, a different 
    # result is obtained for a circular chromosome structure.
    pytest.param(
        [(1,500),(502,1000)], linear_chromosome_struct, (1,1000), id = "lin_gap=1"
    ),
    pytest.param(
        [(1,500),(502,1000)], circular_chromosome_struct, (502,500), id = "circ_gap=1"
    ),
    # a _real_ edge case: all gaps between tuples in the loc_ranges are the 
    # same size (100 seq positions), including the wrap_gap. Since no gaps are
    # larger than wrap_gap, the logic for circular chromosome structure returns
    # the same result as linear.
    pytest.param(
        [(100,199),(300,399),(500,599),(700,799),(900,999)], linear_chromosome_struct, (100,999), id = "lin_equiv_gaps"
    ),
    pytest.param(
        [(100,199),(300,399),(500,599),(700,799),(900,999)], circular_chromosome_struct, (100,999), id = "circ_equiv_gaps"
    ),
    # a _real_ edge case: all gaps between tuples in the loc_ranges are the 
    # same size (100 seq positions), while the wrap_gap size is 99. Since the
    # first gap is larger than wrap_gap, the logic for circular chromosome 
    # structure returns that gap as the separator of start and end positions
    pytest.param(
        [(99,199),(300,399),(500,599),(700,799),(900,999)], linear_chromosome_struct, (99,999), id = "lin_equiv_gaps"
    ),
    pytest.param(
        [(99,199),(300,399),(500,599),(700,799),(900,999)], circular_chromosome_struct, (300,199), id = "circ_equiv_gaps"
    ),
]

@pytest.mark.parametrize("loc_ranges, chr_struct, expected", test_data)
def test_loc_parsing(loc_ranges, chr_struct, expected):
    """ 
    Testing wrapper func for the location parsing function, given a linear 
    chromosome structure type.
    """
    assert parse_embl.process_location_ranges(
        loc_ranges, chr_struct, chromosome_length) == expected


