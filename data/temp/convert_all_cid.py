from all_cid import cid as cid_a
from cid import cid as cid_b

res = {}
for x in cid_a:
    awal = x
    for y in cid_b:
        if(x==y[0]):
            akhir = y[2]
            res[awal] = akhir
            break

print(res)

