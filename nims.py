from math import exp

from drugcipher import Drugcipher
from ppi import PPI

class Nims:
    def __init__(self, compounds, ppi):
        self.ip = self._count_ip()
        self.compounds = compounds
        self.ppi = ppi
        self.proteins = ppi.get_all_proteins()

    def topology_score(self, compA, compB):
        targetA = compA.compound_target
        targetB = compB.compound_target

        tA = ipA = tB = ipB = 0
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
                    if dist > -1:
                        all_dist.append(dist)
                min_dist = min(all_dist)

        for target in targetB:
            if target in targetA:
                min_dist = 0
            else:
                all_dist = []
                for x in targetA:
                    try:
                        dist = self.ppi.prot_dist[(target, x)]
                    except:
                        dist = self.ppi.prot_dist[(x, target)]
                    if dist > -1:
                        all_dist.append(dist)
                min_dist = min(all_dist)
            print(min_dist)


    def agent_score(self, compA, compB):
        pass

    def ip(self, protein=None):
        if protein:
            return self.ip[protein]
        return self.ip

    def _count_ip(self):
        pass

    def synergy_score(self, compA, compB):
        TS = self.topology_score(compA, compB)
        AS = self.agent_score(compA, compB)
        return TS * AS


    # def run(self):
    #     for c in self.compounds:
    #         for o in self.compounds:
    #             if c != o and
    #             (c, o) not in self.all_ss.keys() and
    #             (o, c) not in self.all_ss.keys():
    #                 ss = synergy_score(c, o)
    #                 self.all_ss[(c, o)] = ss

