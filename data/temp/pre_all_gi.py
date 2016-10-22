import json
from all_GI import GI as g
from pre_protein import diabetes as d
from collections import defaultdict

res = defaultdict(list)
total = defaultdict(list)
just_gi = []
for v, k in g: 
    res[v].append(k)

for key, value in res.items():
    if 'More' not in value:
        total[key].append(value)
        just_gi.extend(value)
# print(just_gi)
print(len(res))
# json.dump(res, open("protein_gi.py",'w'))