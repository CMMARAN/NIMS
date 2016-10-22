import json
from protein_gi import protein_gi
from total_prot_481 import total

for t in total:
    t[3] = protein_gi[t[2]]
json.dump(total, open("final_protein_dm.py",'w'))