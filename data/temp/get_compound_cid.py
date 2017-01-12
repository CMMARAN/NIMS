from data.plant_compound_original import plant as plants
from data.cid import cid

res ={}
for plant in plants:
    name = plant[0]
    comps = []
    for i in plant[1]:
        comp = i[0]
        comps.append(comp)
    res[name] = comps
print(res)