from pre_protein import diabetes as d
from total_prot_481 import total as t

prot_d = []
prot_t = []
for i in d:
    if i[1] != '-1':
        prot_d.append(i[1].lower())

for j in t:
    prot_t.append(j[1].lower())

all_entry = []
for i in t:
    all_entry.append(i[2])
print(all_entry)
print(len(all_entry))
