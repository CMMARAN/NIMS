from data.final_prot_dm import protein as prot

class Protein:
    def __init__(self, name, phenotype=None, 
                 gene_name=None, smiles=None):
        self.prot = self._find_prot(name)
        self.name = name
        self.gene_name = self.get_gene_name()
        self.phenotype = self.get_phenotype()
        self.mim = self.get_mim()
        self.gi = self.get_gi()

    def _find_prot(self, name):
        return prot[name]

    def get_name(self):
        return self.prot

    def get_phenotype(self):
        return self.prot[0]

    def get_gene_name(self):
        return None

    def get_mim(self):
        return self.prot[1]

    def get_gi(self):
        return self.prot[2]
