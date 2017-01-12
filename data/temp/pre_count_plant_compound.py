from data.plant_compound_original import plant

count = []
for compound in plant:
    plant_name = compound[0]
    sum_compound = len(compound[1])
    count.append([plant_name, sum_compound])

result = open('data/sum_plant_compound.py','w')
result.write('sum_plant_compound = \n')
result.write(str(count))
