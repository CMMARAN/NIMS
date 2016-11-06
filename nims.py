import numpy as np

from collections import OrderedDict
from math import exp, inf
from sklearn.decomposition import IncrementalPCA
from drugcipher import Drugcipher
from ppi import PPI

class Nims:
    def __init__(self, compounds, ppi):
        self.ppi = ppi
        self.betweenness = ppi.betweenness()
        self.closeness = ppi.closeness()
        self.degree = ppi.degree()
        self.ip = self._count_ip()
        self.compounds = compounds
        self.proteins = ppi.get_all_proteins()


    def _count_TS(self, targetA, targetB):
        tA = ipA = 0
        ip = self.ip

        for target in targetA:
            if target in targetB:
                min_dist = 0
            else:
                all_dist = []
                for x in targetB:
                    try:
                        dist = self.ppi.prot_dist[(target, x)]
                    except:
                        dist = self.ppi.prot_dist[(x, target)]
                    all_dist.append(dist)
                min_dist = min(all_dist)
            tA += exp(-min_dist)*ip[target]
            ipA += ip[target]
        try:
            return tA / ipA
        except:
            return 0


    def topology_score(self, compA, compB):
        targetA = compA.compound_target
        targetB = compB.compound_target

        A = self._count_TS(targetA, targetB)
        B = self._count_TS(targetB, targetA)

        return (A+B) / 2


    def ip(self, protein=None):
        if protein:
            return self.ip[protein]
        return self.ip


    def _count_ip(self):
        ingredient = {}
        b = self.betweenness
        c = self.closeness
        d = self.degree 
        for key, val in b.items():
            ingredient[key] = [val, c[key], d[key]]
        od = OrderedDict(sorted(ingredient.items()))
        all_value = []
        for key, val in od.items():
            all_value.append(val)
        X = np.asarray(all_value)

        ipca = IncrementalPCA(n_components=1, batch_size=3)
        ipca.fit(X)
        res = ipca.transform(X)

        index = 0
        for key in od:
            od[key] = res[index][0] if res[index][0] > 0 else 0
            index += 1
        return od


    def all_TS(self):
        final = {}
        for compA in self.compounds:
            for compB in self.compounds:
                if (compA != compB and
                    (compA.c_id, compB.c_id) not in final.keys()
                    and (compB.c_id, compA.c_id) not in final.keys()):
                    res = self.topology_score(compA, compB)
                    final[(compA.c_id, compB.c_id)] = res
        return final


    def agent_score(self, compA, compB):
        pass


    def synergy_score(self, compA, compB):
        TS = self.topology_score(compA, compB)
        AS = self.agent_score(compA, compB)
        return TS * AS
