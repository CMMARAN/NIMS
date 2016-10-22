import json
from all_cid import cid
from just_gi import just_gi
from senyawa_protein_target import protein


all_target = []
for p in protein:
    if p[0] in cid and p[1] != []:
        all_target.append(p)

target_dm = []
for a in all_target:
    for x in a[1]:
        if x in just_gi:
            target_dm.append(a)
            break
print(target_dm)
# json.dump(target_dm, open("all_target_dm.py",'w'))