
import pytest

import parse_embl

def test_id_line_regex():
    """ Testing wrapper func for the ID line regex pattern. """
    line = "ID   CP002679; SV 1; circular; genomic DNA; STD; PRO; 1038839 BP."
    assert parse_embl.ENA_ID_PATTERN.findall(line) == [("CP002679","circular")]

    line = "ID   BFMR01000110; SV 1; linear; genomic DNA; STD; PRO; 11440 BP."
    assert parse_embl.ENA_ID_PATTERN.findall(line) == [("BFMR01000110","linear")]

    line = "ID   HC710378; SV 1; XXX; protein; PRT; PRO; 409 BP."
    assert parse_embl.ENA_ID_PATTERN.findall(line) == [("HC710378","XXX")]

    line = "FT   source          1..478325"
    assert parse_embl.ENA_ID_PATTERN.findall(line) == []

def test_FT_block_line_regex():
    """ Testing wrapper func for the Feature Key line pattern. """
    lines = """ID   ABZA01000001; SV 1; linear; genomic DNA; WGS; PRO; 478325 BP.
XX
FT   source          1..478325
FT                   /organism="Wolbachia endosymbiont of Culex quinquefasciatus
FT                   JHB"
FT                   /db_xref="taxon:569881"
FT   gene            <1..1701
FT                   /locus_tag="C1A_288"
FT   CDS             <1..1701"""
    ground_truth = [False, False, True, False, False, False, True, False, True]
    matches = [True if parse_embl.FT_START_PATTERN.findall(line) else False 
               for line in lines.split('\n')]
    assert matches == ground_truth

def test_XREF_lines_regex():
    """ Testing wrapper func for the XREF line patterns. """
    lines = """FT   CDS             <1..1701
FT                   /db_xref="InterPro:IPR023614"
FT                   /db_xref="UniProtKB/TrEMBL:B6Y618"
FT                   /protein_id="EEB56106.1"
FT   CDS             complement(1822..1956)
FT                   /locus_tag="C1A_289"
FT                   /db_xref="UniProtKB/TrEMBL:B6Y619"
FT                   /protein_id="EEB56107.1"
FT                   /translation="MLKYNVSDDDGKMDPSVKHWDDTIYYANCHNFRTAVTGMTLLIV" """
    ground_truth = [False, False, True, True, False, False, True, True, False]
    matches = [True if parse_embl.XREF_SEARCH_PATTERN.findall(line) else False for line in lines.split('\n')]
    assert matches == ground_truth

def test_location_lines_regex():
    """ Testing wrapper func for the Location line patterns. """
    lines = """FT   CDS             J00194.1:100..202
FT   CDS             467        # we ignore these
FT   CDS             340..565
FT   CDS             <345..500
FT   CDS             <1..888
FT   CDS             1..>888
FT   CDS             102.110    # we ignore these
FT   CDS             123^124    # we ignore these
FT   CDS             join(12..78,134..202)
FT   CDS             complement(34..126)
FT   CDS             complement(join(2691..4571,4918..5163))
FT   CDS             join(complement(4918..5163),complement(2691..4571))
FT   CDS             join(1..100,J00194.1:100..202) """
    ground_truth = [
        [('100','202')],
        [],
        [('340','565')],
        [('345','500')],
        [('1','888')],
        [('1','888')],
        [],
        [],
        [('12','78'), ('134','202')],
        [('34','126')],
        [('2691','4571'),('4918','5163')],
        [('4918','5163'),('2691','4571')],
        [('1','100'),('100','202')]
    ]
    locs = [parse_embl.CDS_LOC_PATTERN.findall(line) for line in lines.split('\n')]
    print(locs)
    assert locs == ground_truth


