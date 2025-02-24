
import pickle
import sys

tab_in = sys.argv[1]
pkl_out= sys.argv[2]

id_dict = {}
with open(tab_in,'r') as tab:
    for line in tab:
        line_list = line.strip().split()
        id_ = line_list[0]
        uniprot_id = line_list[1]
        rest = line_list[2:]
        subdict = id_dict.get(id_,{})
        subdict[uniprot_id] = rest
        id_dict[id_] = subdict

with open(pkl_out,'wb') as tab_pkl:
    pickle.dump(["CONTENTS: [0] this string and [1] dict of ENA Ids mapping to dicts of UniProtIds mapping to lists of the tab fields",id_dict],tab_pkl)

