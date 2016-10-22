from all_cid import cid
from protein_list import prot

res = {}
x = 0
for c in cid:
    target = []
    for p in prot:
        if c in p[1]:
            target.append(p[0])
    if(len(target)!=0):
        x+=1
        print(x)
    res[c] = target
