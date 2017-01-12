from flask import Flask, render_template
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
app = Flask(__name__)


@app.route("/")
def index():
    return render_template("main.html", data=compounds)


if __name__ == "__main__":
    app.run(port=5000, host="localhost")
