from gi_sisa import gi_sisa as g
from collections import defaultdict

res = defaultdict(list)
for v, k in g: 
    res[v].append(k)
print(res)