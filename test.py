from compound import Compound
from data.ppi_data import ppi
from data.comp_cid import cid
from drugcipher import Drugcipher
from nims import Nims
from ppi import PPI

compounds = []
for c in cid:
    senyawa = Compound(c)
    compounds.append(senyawa)

g = PPI("diabetes mellitus", ppi)

drug = Drugcipher(ppi=g, compounds=compounds)
drug.run()

net = Nims(ppi=g, compounds=compounds)
net.all_TS()
