import csv
import json
from final_prot_dm import protein as p
from collections import defaultdict

final_dict = defaultdict(list)

with open("test.tsv", 'r') as f:
    gene_pheno = csv.reader(f, delimiter='\t')
    for row in gene_pheno:
        if row[1] in p.keys():
            final_dict[row[1]].append(row[3])

print(final_dict)
print(len(final_dict))
json.dump(final_dict, open("final_gene_phenotype.py",'w'))
