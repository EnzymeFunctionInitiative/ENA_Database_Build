
import gzip

from Bio import SeqIO

###############################################################################
# Parsing an EMBL flat file
###############################################################################

def process_file(
        file_path: str,
        database_connection,
        output_file: str):
    """
    """

    process_results = {}
    # open and read the gzipped file
    with gzip.open(file_path, 'rt') as handle:
        # use biopython to parse the embl flat file
        # loop over the ENA accession Id'd records in the file
        for record in SeqIO.parse(handle, "embl"):
            # grab the record's ENA accession name
            ID = record.name
            # consider the chromosome/genetic topology type
            if record.annotations['topology'] in ["linear","circular"]:
                CHR = 1 if record.annotations['topology'] == "linear" else 0
            else:
                print("!!! Unknown chromosome type " 
                        + f"{record.annotations['topology']} observed in "
                        + f"{file_path}: {record.name}")
                ID = ""
            
            # start a fresh count for this ID
            count = 0
            # loop over features associated with the record
            for feature in record.features:
                # only consider coding sequence features
                if feature.type != "CDS":
                    continue
                
                count += 1
                # grab the start and end locations for the feature
                # add one to both because biopython uses pythonic zero indexing
                START = int(feature.location.start) + 1
                END   = int(feature.location.end)
                # grab the stand direction for the feature
                # feature.location.strand == -1 if "complement" is in line 
                # else == 1. Our current expected format is 0 if complement
                # else 1.
                DIR = 0 if feature.location.strand == -1 else 1
                # grab the protein_ids
                if feature.qualifiers['protein_id']:
                    proteinIds = [protId.split('.')[0] for protId in feature.qualifiers['protein_id']]
                # grab the uniprot_ids
                if feature.qualifiers['db_xref']:
                    uniprotIds = [xref.split(':')[1] for xref in feature.qualifiers['db_xref'] if 'UniProtKB' in xref]

                # do database reverseLookup
                rev_uniprot_ids, no_match = database_connection.reverse_lookup(proteinIds)
                # check whether the rev_uniprot_ids set is filled with 
                # uniprotIds; if it is, use those as the IDs to IO to file
                if rev_uniprot_ids:
                    uniprotIds = rev_uniprot_ids
            
                # loop over uniprotIds and append to file. If uniprotIds is 
                # empty, the file is not written to.
                for id_ in uniprotIds:
                    with open(output_file,"a") as out_tab:
                        # format: {ENA_ID}\t{UniprotId}\t{seqCount}\t{topologyInt}\t{directionInt}\t{START}\t{END}\n
                        out_tab.write(f"{ID}\t{id_}\t{count}\t{CHR}\t{DIR}\t{START}\t{END}\n")
    
    return output_file
