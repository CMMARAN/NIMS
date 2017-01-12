import csv
from data.cid import cid

res = []
for m in cid:
    if m[2] == "-1":
        continue
    res.append(int(m[2]))

result = open('data/cid_only.py','w')
result.write('cid = \n')
result.write(str(res))
