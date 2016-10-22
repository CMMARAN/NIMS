import json
from final_compound_dm import compound
from final_protein_dm import protein

res = {}
for key, value in compound.items():
    target = []
    for x in value[6]:
        for p in protein:
            if x in p[3]:
                target.append(p[1])
    value[6] = target
