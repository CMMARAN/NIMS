import networkx as nx
import numpy as np
from indigo.indigo import *
from compound import Compound
from protein import Protein
from data.ppi_data import ppi
from data.prot_name import prot_name as prot
from data.comp_cid import cid
from drugcipher import Drugcipher
from ppi import PPI
from math import exp

def count_prot_dist_all(g, prots):
    res = {}
    for p in prots:
        for r in prots:
            if p != r and ((p,r) not in res.keys() and (r,p) not in res.keys()):
                distance = float(g.shortest_path(p, r))
                res[(p,r)] = distance
    return res

def cs_all(compound, compounds=None):
    row = []
    for compB in compounds:
        cs = count_cs(compound, compB)
        row.append(cs)
    return row

def prot_comp_closeness_all(protein, compounds):
    dist = []
    for compound in compounds:
        compound_targets = compound.get_compound_target()
        total_distance = prot_comp_closeness(protein, compound_targets)
        dist.append(total_distance)
    return dist

def prot_comp_closeness(protein, compound_targets):
    total_distance = 0
    for compound_target in compound_targets:
        distance = float(g.shortest_path(protein, compound_target))
        if distance != -1.2:
            closeness = exp((distance*(-1))**2)
        else:
            closeness = 0
        total_distance += closeness
    return total_distance

g = PPI("diabetes mellitus", ppi)
s = Compound("C00001432")

prots = g.get_all_proteins()
res = count_prot_dist_all(g, prots)
print(res)
print(len(res))

# target = s.compound_target
# print(target)
# prot = Protein("IFI27")

# compounds = []
# all_target = []
# for c in cid:
#     senyawa = Compound(c)
#     compounds.append(senyawa)
#     for x in senyawa.compound_target:
#         if x not in all_target:
#             all_target.append(x)
# print(all_target)
# print(len(all_target))
# key = [x for x in g.shortest_path("IFI27").keys()]
# z = []
# for x in key:
#     if x not in all_target:
#         z.append(x)
# print(z)
# print(key)
# print(len(z))
# print(len(key))

# compounds = []
# for c in cid:
#     senyawa = Compound(c)
#     compounds.append(senyawa)
# cs = cs_all(s, compounds)

# print(cs)
# prot_comp = prot_comp_closeness_all("IFI27", compounds)
# print(prot_comp)
# if all(v == 0 for v in prot_comp):
#     print('indeed they are')
# concordance = np.cov(cs,prot_comp) / (np.std(cs) * np.std(prot_comp))

# print(concordance)



