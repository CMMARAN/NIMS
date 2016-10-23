from final_protein_dm import protein
import json

res = {}
for p in protein:
    value = []
    value.append(p[0])
    value.append(p[2])
    value.append(p[3])
    res[p[1]] = value

json.dump(res, open('final_prot_dm.py', 'w'))