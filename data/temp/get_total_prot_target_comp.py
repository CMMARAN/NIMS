from data.final_comp_dm_set import compounds as cp
from data.plant_compound_original import plant as pl

for plant, compounds in pl:
    prot = []
    for comp in compounds:
        if comp[0] in cp.keys():
            data = cp[comp[0]]
            prot = prot + data[6]
            total = set(prot)
    print(plant, len(total))