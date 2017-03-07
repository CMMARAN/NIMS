import numpy as np
import pickle
from collections import OrderedDict
from math import exp
from sklearn.decomposition import IncrementalPCA

from data.data_as import pheno_sim
from data.final_gene_phenotype import gene_phenotype as pheno_db
from hpo_similarity.hpo_similarity.ontology import Ontology
from hpo_similarity.hpo_similarity.check_proband_terms import check_terms_in_graph
from hpo_similarity.hpo_similarity.test_similarity import test_similarity


class Nims:

    def __init__(self, compounds, ppi):
        self.ppi = ppi
        self.betweenness = ppi.betweenness()
        self.closeness = ppi.closeness()
        self.degree = ppi.degree()
        self.ip = self._count_ip()
        self.compounds = compounds
        self.proteins = ppi.get_all_proteins()
        self.hpo_ontology = Ontology()
        self.pheno_graph = None
        self.alt = self.hpo_ontology.get_alt_ids()
        self.obsolete = self.hpo_ontology.get_obsolete_ids()
        self.pheno_db = pheno_db
        self.AS = self.load_AS()
        self.SS = self.load_SS()
        self.TS = self.load_TS()
        self.pheno_sim = pheno_sim

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
            tA += exp(-min_dist) * ip[target]
            ipA += ip[target]
        try:
            return tA / ipA
        except:
            return 0

    def topology_score(self, compA, compB):
        if ((compA.metabolite, compB.metabolite) in self.TS.keys() or
                (compB.metabolite, compA.metabolite) in self.TS.keys()):
            try:
                TS = self.TS[(compA.metabolite, compB.metabolite)]
            except:
                TS = self.TS[(compB.metabolite, compA.metabolite)]
            return TS

        targetA = compA.compound_target
        targetB = compB.compound_target

        A = self._count_TS(targetA, targetB)
        B = self._count_TS(targetB, targetA)

        return (A + B) / 2

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
            ingredient[key].append(od[key])
            index += 1
        return od

    def all_TS(self):
        for compA in self.compounds:
            for compB in self.compounds:
                if (((compA.metabolite, compB.metabolite) not in self.TS.keys() or
                    (compB.metabolite, compA.metabolite) not in self.TS.keys())):
                    res = self.topology_score(compA, compB)
                    self.TS[(compA.metabolite, compB.metabolite)] = res
                    self.TS[(compB.metabolite, compA.metabolite)] = res

        pickle.dump(self.TS, open("data/TS1.pkl", "wb"))
        return self.TS

    def generate_HPO(self):
        obs = list(self.obsolete)
        self.pheno_graph = self.hpo_ontology.get_graph()
        for prots in self.pheno_db:
            terms = self.pheno_db[prots]
            terms = [self.alt[term]
                     if term in self.alt else term for term in terms]
            terms1 = [term for term in terms if term not in obs]
            self.pheno_db[prots] = terms1
        self.pheno_graph.tally_hpo_terms(self.pheno_db)
        check_terms_in_graph(self.pheno_graph, self.pheno_db)

    def generate_phenotype_similarity(self, protA, protB):
        if protA not in self.pheno_db.keys() or protB not in self.pheno_db.keys():
            return 0

        proteins = [protA, protB]
        PS = test_similarity(
            self.pheno_graph, self.pheno_db, proteins, "resnik")
        return PS if PS > 0 else 0

    def agent_score(self, compA, compB):
        if ((compA.metabolite, compB.metabolite) in self.AS.keys() or
                (compB.metabolite, compA.metabolite) in self.AS.keys()):
            try:
                AS = self.AS[(compA.metabolite, compB.metabolite)]
            except:
                AS = self.AS[(compB.metabolite, compA.metabolite)]
            return AS

        targetA = compA.compound_target
        targetB = compB.compound_target
        count = 0
        total_PS = 0
        for x in targetA:
            for y in targetB:
                if ((x, y) not in self.pheno_sim.keys() and
                        (y, x) not in self.pheno_sim.keys()):
                    PS = self.generate_phenotype_similarity(x, y)
                    self.pheno_sim[(x, y)] = PS
                    self.pheno_sim[(y, x)] = PS
                else:
                    try:
                        PS = self.pheno_sim[(x, y)]
                    except:
                        PS = self.pheno_sim[(y, x)]
                total_PS += PS
                count += 1

        return total_PS / count if count > 0 else 0

    def all_AS(self):
        self.generate_HPO()
        for compA in self.compounds:
            for compB in self.compounds:
                if (((compA.metabolite, compB.metabolite) not in self.AS.keys() or
                        (compA.metabolite, compB.metabolite) not in self.AS.keys())):
                    res = self.agent_score(compA, compB)
                    self.AS[(compA.metabolite, compB.metabolite)] = res
                    self.AS[(compB.metabolite, compA.metabolite)] = res

        pickle.dump(self.AS, open("data/AS.pkl", "wb"))
        return self.AS

    def synergy_score(self, compA, compB):
        if ((compA.metabolite, compB.metabolite) in self.SS.keys() or
                (compB.metabolite, compA.metabolite) in self.SS.keys()):
            try:
                SS = self.SS[(compA.metabolite, compB.metabolite)]
            except:
                SS = self.SS[(compB.metabolite, compA.metabolite)]
            return SS

        TS = self.topology_score(compA, compB)
        AS = self.agent_score(compA, compB)
        return TS * AS

    def all_SS(self):
        for compA in self.compounds:
            for compB in self.compounds:
                if (((compA.metabolite, compB.metabolite) not in self.SS.keys() or
                        (compA.metabolite, compB.metabolite) not in self.SS.keys())):
                    res = self.synergy_score(compA, compB)
                    self.SS[(compA.metabolite, compB.metabolite)] = res
                    self.SS[(compB.metabolite, compA.metabolite)] = res

        pickle.dump(self.SS, open("data/SS.pkl", "wb"))
        return self.SS

    def load_TS(self):
        with open("data/TS.pkl", "rb") as f:
                self.TS = pickle.load(f)
        return self.TS

    def load_AS(self):
        with open("data/AS.pkl", "rb") as f:
                self.AS = pickle.load(f)
        return self.AS

    def load_SS(self):
        with open("data/SS.pkl", "rb") as f:
                self.SS = pickle.load(f)
        return self.SS
