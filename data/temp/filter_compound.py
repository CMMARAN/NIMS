import json
from dm_comps import comps
from final_target_dm import target

res = {}
for t in target:
    comps[t[0]].append(t[1])
    res[t[0]] = comps[t[0]]
print(res)
json.dump(res, open("final_compound_dm.py",'w'))
