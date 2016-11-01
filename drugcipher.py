import numpy as np
import operator
from time import time, mktime
from ppi import PPI
from compound import Compound
from indigo.indigo import *
from math import exp
from collections import OrderedDict


class Drugcipher:
    def __init__(self, ppi, compounds):
        self.ppi = ppi
        self.compounds = compounds
        self.proteins = ppi.get_all_proteins()
        self.prot_dist = ppi.prot_dist
        self.chem_sim = {}


    def count_cs(self, compA, compB):
        c = self.chem_sim
        if (compA.c_id, compB.c_id) not in c.keys() and (compB.c_id, compA.c_id) not in c.keys():
            indigo = Indigo()
            m1 = indigo.loadMolecule(compA.smiles)
            m2 = indigo.loadMolecule(compB.smiles)
            m1.aromatize()
            m2.aromatize()
            fp1 = m1.fingerprint("sim")
            fp2 = m2.fingerprint("sim")
            cs = indigo.similarity(fp1, fp2, "tanimoto")
            self.chem_sim[(compA.c_id, compB.c_id)] = cs
        else:
            try:
                cs = c[(compA.c_id, compB.c_id)]
            except:
                cs = c[(compB.c_id, compA.c_id)]

        return cs


    def cs_all(self, compound, compounds=None):
        row = []
        for compB in compounds:
            cs = self.count_cs(compound, compB)
            row.append(cs)
        return row


    def prot_comp_closeness(self, protein, compound_targets):
        total_distance = 0
        key = self.prot_dist.keys()
        for c in compound_targets:
            if protein == c:
                distance = 0
            else:
                if (protein, c) not in key and (c, protein) not in key:
                    distance = float(self.ppi.shortest_path(protein, c)) 
                    self.prot_dist[(protein, c)] = distance
                else:
                    try:
                        distance = self.prot_dist[(protein, c)]
                    except:
                        distance = self.prot_dist[(c, protein)]
            if distance != -1.2:
                closeness = exp((distance*(-1))**2)
            else:
                closeness = 0
            total_distance += closeness
        return total_distance


    def prot_comp_closeness_all(self, proteins, compounds):
        res = {}
        for protein in proteins:
            dist = []
            for compound in compounds:
                compound_targets = compound.get_compound_target()
                total_distance = self.prot_comp_closeness(protein, compound_targets)
                dist.append(total_distance)
            res[protein] = dist

        return res

    def cov(self, a, b):
        if len(a) != len(b):
            return
        a_mean = np.mean(a)
        b_mean = np.mean(b)
        sum = 0
        for i in range(0, len(a)):
            sum += ((a[i] - a_mean) * (b[i] - b_mean))
        return sum/(len(a)-1)


    def get_proteins_rank(self, compound, limit):
        rank = {}
        cs = self.cs_all(compound=compound, compounds=self.compounds)
        prot_comp = self.prot_comp_closeness_all(self.proteins, self.compounds)
        for prot, distance in prot_comp.items():
            if all(v == 0 for v in distance):
                concordance = 0
            else:
                concordance = self.cov(cs, distance) / (np.std(cs) * np.std(distance))
            rank[prot] = concordance
        rank_ordered = sorted(rank.items(), key=operator.itemgetter(1), reverse=True)
        prots = []
        i = 0
        for x in rank_ordered:
            if i < limit:
                prots.append(x[0])
            else:
                break
            i+=1
        print("Proteins to be added to", compound.c_id,":", prots)
        return prots


    def extend_compound_target(self, compound):
        print("Starting extend compound: ", compound.c_id)
        compound_targets = compound.get_compound_target()
        print("Compound target: ", compound_targets)
        proteins = self.get_proteins_rank(compound, 5)
        for p in proteins:
            if p not in compound_targets:
                compound.set_compound_target(p)
        compound_targets = compound.get_compound_target()
        print("Compound target now", compound_targets)


    def run(self):
        awal = time()
        print("Total protein distance awal:", len(self.prot_dist))
        for compound in self.compounds:
            self.extend_compound_target(compound)
        akhir = time()
        gap = akhir - awal
        print("Total chemical similiarty:", len(self.chem_sim))
        print("Total protein distance akhir:", len(self.prot_dist))
        print("Script done in", gap, "seconds")
