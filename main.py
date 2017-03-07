import pickle
import operator
from collections import defaultdict
from math import ceil
from flask import Flask, render_template, request, url_for
from compound import Compound
from data.comp_cid import cid

compounds = []
for c in cid:
    senyawa = Compound(c)
    compounds.append(senyawa)

app = Flask(__name__)

with open("data/cs.pkl", "rb") as f:
    CS = pickle.load(f)

with open("data/TS.pkl", "rb") as f:
    TS = pickle.load(f)

with open("data/AS.pkl", "rb") as f:
    AS = pickle.load(f)

with open("data/SS.pkl", "rb") as f:
    SS = pickle.load(f)
    sorted_SS = sorted(SS.items(), key=operator.itemgetter(1), reverse=True)
    unique_SS = {}
    rank = 0
    for idx, val in enumerate(sorted_SS):
        if idx % 2 == 0:
            rank += 1
            unique_SS[val[0]] = (round(val[1], 4), rank)


def _get_navigation_url():
    urls = [
        ("NIMS", url_for("index")),
        ("All Ranking", url_for("synergy")),
        ("About", url_for("about"))
    ]
    return urls


def _get_score(compA, compB):
    if (compA, compB) not in CS.keys() and (compA, compB) not in CS.keys():
        chem = "Not found"
    else:
        try:
            chem = "%.4f" % CS[(compA, compB)]
        except:
            chem = "%.4f" % CS[(compB, compA)]

    if (compA, compB) not in TS.keys() and (compA, compB) not in TS.keys():
        topology = "Not found"
    else:
        try:
            topology = "%.4f" % TS[(compA, compB)]
        except:
            topology = "%.4f" % TS[(compB, compA)]

    if (compA, compB) not in TS.keys() and (compA, compB) not in AS.keys():
        agent = "Not found"
    else:
        try:
            agent = "%.4f" % AS[(compA, compB)]
        except:
            agent = "%.4f" % AS[(compB, compA)]

    if (compA, compB) not in TS.keys() and (compA, compB) not in SS.keys():
        synergy = "Not found"
        rank = "Not found"
    else:
        try:
            synergy = "%.4f" % SS[(compA, compB)]
            rank = sorted_SS.index(((compA, compB), SS[(compA, compB)]))
            rank = "{} of {}".format(ceil((rank + 1) / 2), int(len(SS) / 2))
        except:
            synergy = "%.4f" % SS[(compB, compA)]
            rank = sorted_SS.index(((compB, compA), SS[(compB, compA)]))
            rank = "{} of {}".format(ceil((rank + 1) / 2), int(len(SS) / 2))

    targets = defaultdict(list)
    for c in compounds:
        if c.metabolite == compA:
            targets[compA] = c.compound_target
        if c.metabolite == compB:
            targets[compB] = c.compound_target

    data = {
        "topology": topology,
        "agent": agent,
        "synergy": synergy,
        "chem": chem,
        "rank": rank,
        "targets": targets
    }
    return data


@app.route("/", methods=['POST', 'GET'])
def index():
    nav = _get_navigation_url()
    comp_data = {
        "topology": "",
        "agent": "",
        "synergy": "",
        "chem": "",
        "rank": "",
        "targets": {"": "", " ": ""}
    }

    if(request.method == 'POST'):
        compA = request.form['compA']
        compB = request.form['compB']
        comp_data = _get_score(compA, compB)
    return render_template("main.html",
                           data=compounds,
                           comp_data=comp_data,
                           navigation_url=nav)


@app.route("/about", methods=['POST', 'GET'])
def about():
    nav = _get_navigation_url()
    return render_template("about.html",
                           navigation_url=nav)


@app.route("/synergy", methods=['POST', 'GET'])
def synergy():
    nav = _get_navigation_url()
    return render_template("synergy.html",
                           data=unique_SS,
                           navigation_url=nav)


if __name__ == "__main__":
    app.run(port=5000, host="localhost")
