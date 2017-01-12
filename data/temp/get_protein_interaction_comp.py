from data.final_comp_dm_set import compounds as cp
from data.plant_compound_original import plant as pl

for plant, compounds in pl:
    total = 0
    for comp in compounds:
        if comp[0] in cp.keys():
            data = cp[comp[0]]
            total += len(data[6])
    print(plant, total)