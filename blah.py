from data.final_gene_phenotype import gene_phenotype as g

count = 0
for x, y in g.items():
    count += len(y)
print(count)