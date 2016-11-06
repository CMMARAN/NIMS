from data.final_comp_dm_set import compounds as comps

class Compound:
    def __init__(self, c_id, metabolite=None, cas_id=None,
                 molecul_formula=None, mass=None, plants=None,
                 smiles=None, compound_target=None):
        self.c_id = c_id
        self.comp = self._find_compound(c_id)
        self.cas_id = self.get_cas_id()
        self.metabolite = self.get_metabolite()
        self.molecul_formula = self.get_molecul_formula()
        self.mass = self.get_mass()
        self.plants = self.get_plants()
        self.smiles = self.get_smiles()
        self.compound_target = self.get_compound_target()

    def _find_compound(self, c_id):
        return comps[c_id]
            

    def get_c_id():
        return self.c_id

    def get_cas_id(self):
        return self.comp[0]

    def get_metabolite(self):
        return self.comp[1]

    def get_molecul_formula(self):
        return self.comp[2]

    def get_mass(self):
        return self.comp[3]

    def get_plants(self):
        return self.comp[4]

    def get_smiles(self):
        return self.comp[5]

    def get_compound_target(self):
        return self.comp[6]

    def set_compound_target(self, protein):
        if protein not in self.comp[6]:
            self.compound_target.append(protein)
        return
