import json
from all_target_dm import all_target_dm
from just_gi import just_gi

final = []
i = 0
for a in all_target_dm:
    tmp = len(a[1])
    res = []
    for x in a[1]:
        if x in just_gi:
            # print(x)
            res.append(x)
    a[1] = res
    final.append(a)
    if len(res) != 0:
        i+=1
    print(a[0], tmp, len(a[1]))
print(final)
json.dump(final, open("final_target_dm.py",'w'))
