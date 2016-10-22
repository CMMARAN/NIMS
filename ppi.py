import networkx as nx

class PPI:
    def __init__(self, disease, graph):
        self.disease = disease
        self.graph = graph
        self.network = self._network()


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
        return nx.shortest_path_length(self.network, source=proteinA, target=proteinB)

    def get_all_proteins(self):
        n = self.network
        return n.nodes()
