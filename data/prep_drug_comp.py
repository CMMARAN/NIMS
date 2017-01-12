from final_comp_dm_set import compounds
from drugcipher_result import drug_res

for key, res in drug_res.items():
    compounds[key][6] = res

result = open('final_comp_dm_drugcipher.py','w')
result.write('compounds = \n')
result.write(str(compounds))