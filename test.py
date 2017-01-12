from compound import Compound
from data.ppi_data import ppi
from data.comp_cid import cid
from nims import Nims
from ppi import PPI

compounds = []
for c in cid:
    senyawa = Compound(c)
    compounds.append(senyawa)

g = PPI("diabetes mellitus", ppi)


net = Nims(ppi=g, compounds=compounds)
