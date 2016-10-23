import networkx as nx

class PPI:
    def __init__(self, disease, graph):
        self.disease = disease
        self.graph = graph
        self.network = self._network()
        self.prot_dist = self.get_all_shortest_path()


    def name(self):
        return self.disease


    def _network(self):
        g = nx.Graph()
        g.add_edges_from(self.graph)
        return g


    def betweenness(self, protein=None):
        b = nx.betweenness_centrality(self.network)
        if protein:
            protein = protein.upper()
            for key in b:
                if key == protein:
                    return b[protein]
            return ("%s not found in this PPI"%protein)
        return sorted(b.items(), key=lambda x: x[1], reverse=True)


    def closeness(self, protein=None):
        c = nx.closeness_centrality(self.network)
        if protein:
            protein = protein.upper()
            for key in c:
                if key == protein:
                    return c[protein]
            return ("%s not found in this PPI"%protein)
        return sorted(c.items(), key=lambda x: x[1], reverse=True)


    def degree(self, protein=None):
        d = nx.degree_centrality(self.network)
        if protein:
            protein = protein.upper()
            for key in d:
                if key == protein:
                    return d[protein]
            return ("%s not found in this PPI"%protein)
        return sorted(d.items(), key=lambda x: x[1], reverse=True)


    def shortest_path(self, proteinA=None, proteinB=None):
        try:
            shortest = nx.shortest_path_length(self.network, source=proteinA, target=proteinB)
        except:
            shortest = -1.2
        return shortest

    def get_all_proteins(self):
        n = self.network
        return n.nodes()

    def get_all_shortest_path(self):
        res = {}
        prots = self.get_all_proteins()
        for p in prots:
            for r in prots:
                if p != r and ((p,r) not in res.keys() and (r,p) not in res.keys()):
                    distance = float(self.shortest_path(p, r))
                    res[(p,r)] = distance
        return res
