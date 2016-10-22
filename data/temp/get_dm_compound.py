from compound_list import comps
from cid import cid as CID
from all_cid import cid
from smiles1 import smiles

import json

res = {}
for c in cid:
    for comp in comps:
        if c == comp[0]:
            temp = comp[1]
            for x in CID:
                if c == x[0]:
                    if x[2] != '-1':
                        temp.append(smiles[x[2]])
                        res[c] = temp
                        break
json.dump(res, open("blah.py",'w'))