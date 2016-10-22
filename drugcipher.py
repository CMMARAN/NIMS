import numpy as np

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


    @staticmethod
    def count_cs(compA, compB):
        #hitung tanimoto
        indigo = Indigo()
        m1 = indigo.loadMolecule(compA.smiles)
        m2 = indigo.loadMolecule(compB.smiles)
        m1.aromatize()
        m2.aromatize()
        fp1 = m1.fingerprint("sim")
        fp2 = m2.fingerprint("sim")
        cs = indigo.similarity(fp1, fp2, "tanimoto")

        return cs


    @staticmethod
    def cs_all(compound, compounds=None):
        row = []
        for compB in compounds:
            cs = count_cs(compound, compB)
            row.append(cs)
        return CS


    @staticmethod
    def prot_comp_closeness(protein, compound_targets):
        total_distance = 0
        for compound_target in compound_targets:
            distance = ppi.shortest_path(protein, compound_target)
            closeness = exp(distance*(-1))
            total_distance += closeness
        return total_distance


    def prot_comp_closeness_all(proteins, compounds):
        res = {}
        for protein in proteins:
            dist = []
            for compound in compounds:
                compound_targets = compound.get_compound_target()
                total_distance = prot_comp_closeness(protein, compound_targets)
                dist.append(total_distance)
            res[protein] = dist
        return res


    def get_proteins_rank(compound, limit):
        rank = {}
        cs = cs_all(compound=compound, compounds=self.compounds)
        prot_comp = prot_comp_closeness_all(self.proteins, self.compounds)
        for prot, distance in prot_comp:
            concordance = np.cov(cs,distance) / (np.std(cs) * np.std(distance))
            rank[prot] = concordance
        rank_ordered = OrderedDict(sorted(rank.items(), reverse=True))
        prots = []
        i = 0
        for prot, value in rank_ordered:
            if i < limit:
                prots.append(prot)
            else:
                break
            i+=1
        return prots


    def extend_compound_target(compound):
        compound_targets = compound.get_compound_target()
        proteins = self.get_proteins_rank(compound, 100)
        for p in proteins:
            if p not in compound_targets:
                compound.set_compound_target(p)


    def run():
        for compound in self.compounds:
            extend_compound_target(compound)
