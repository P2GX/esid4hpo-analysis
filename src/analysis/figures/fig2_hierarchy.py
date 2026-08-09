#!/usr/bin/env python3
"""
Figure 2 - worked structural examples (before/after hierarchy diagrams).

Each example is a hand-curated set of nodes and edges (you don't show the whole
subtree - you pick the terms that make the change legible) rendered as a
before/after pair with Graphviz. Node kinds are colour-coded:

    root   BIH navy            the branch anchor
    new    BIH blue            term created during the initiative
    rep    grey + amber edge   pre-existing term re-parented into the new scaffold
    exist  grey                pre-existing term, unchanged
    ctx    dashed outline      context (a parent branch shown for orientation)

The node/edge sets below were extracted from the real releases
(v2024-08-13 vs v2026-06-23) - verify against the current edit file before
finalising, e.g. with hp-edit.owl in Protege or:
    grep "HP_5210337" hp-edit.owl | grep SubClassOf

Deps:  pip install graphviz   (+ system graphviz: `brew install graphviz`)
Run:   python fig2_hierarchy.py
Outputs: _work/figures/fig2_<example>_{before,after}.{svg,pdf,png}
"""

import pathlib
from graphviz import Digraph

OUTDIR = pathlib.Path(__file__).resolve().parent / "_work" / "figures"

NAVY="#003754"; BLUE="#4e7e96"; GREY="#e7eaec"; INK="#2b2f33"
AMBER="#c8a870"; EDGE="#9aa3ad"; BORDER="#c4c9ce"


def _node(g, nid, label, kind):
    lab = label.replace(" ", "\n") if len(label) > 16 else label
    styles = {
        "root":  dict(fillcolor=NAVY, fontcolor="white", color=NAVY),
        "new":   dict(fillcolor=BLUE, fontcolor="white", color=BLUE),
        "rep":   dict(fillcolor=GREY, fontcolor=INK, color=AMBER, penwidth="2.4"),
        "exist": dict(fillcolor=GREY, fontcolor=INK, color=BORDER),
        "ctx":   dict(fillcolor="white", fontcolor="#9aa3ad", color=BORDER,
                      style="rounded,dashed"),
    }
    g.node(nid, lab, **styles[kind])


def _graph(name):
    g = Digraph(name)
    g.attr(rankdir="TB", bgcolor="white", splines="polyline")
    g.attr("node", shape="box", style="rounded,filled", fontname="Helvetica",
           fontsize="12", margin="0.14,0.07", penwidth="1.1")
    g.attr("edge", color=EDGE, penwidth="1.2", arrowsize="0.7")
    return g


def render(example, nodes, edges):
    """nodes: {nid:(label,kind)}; edges: [(child,parent,{amber?})]"""
    g = _graph(example)
    for nid, (label, kind) in nodes.items():
        _node(g, nid, label, kind)
    for e in edges:
        child, parent = e[0], e[1]
        amber = len(e) > 2 and e[2]
        g.edge(child, parent, color=AMBER if amber else EDGE,
               penwidth="1.8" if amber else "1.2")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for fmt in ("svg", "pdf", "png"):
        g.render(str(OUTDIR / example), format=fmt, cleanup=True)
    print(f"rendered {example}")


# ============================================================================
# EXAMPLE 1 - granuloma / granulomatosis etiology axis  (verified v2024 vs v2026)
# ============================================================================
GRANULOMA_BEFORE = (
    {"macro": ("Abnormal macrophage morphology", "ctx"),
     "pulm_i": ("Abnormal pulmonary interstitial morphology", "ctx"),
     "gran": ("Granuloma", "ctx"),
     "G": ("Granulomatosis", "exist"),
     "PG": ("Pulmonary granulomatosis", "exist"),
     "NPG": ("Necrotizing pulmonary granulomatosis", "exist"),
     "NNPG": ("Non-necrotizing pulmonary granulomatosis", "exist"),
     "EG": ("Eosinophilic granuloma", "exist")},
    [("G", "macro"), ("PG", "pulm_i"), ("EG", "gran"),
     ("NPG", "PG"), ("NNPG", "PG")],
)
GRANULOMA_AFTER = (
    {"G": ("Granulomatosis", "root"),
     "INF": ("Infectious granulomatosis", "new"),
     "NINF": ("Non-infectious granulomatosis", "new"),
     "CAS": ("Caseating granulomatosis", "new"),
     "NCAS": ("Non-caseating granulomatosis", "new"),
     "GLILD": ("Granulomatous lymphocytic interstitial lung disease", "new"),
     "GLN": ("Granulomatous lymphadenopathy", "new"),
     "GA": ("Granuloma annulare", "new"),
     "NPG": ("Necrotizing pulmonary granulomatosis", "rep"),
     "NNPG": ("Non-necrotizing pulmonary granulomatosis", "rep"),
     "EG": ("Eosinophilic granuloma", "rep")},
    [("INF", "G"), ("NINF", "G"), ("GLN", "G"),
     ("CAS", "INF"), ("NCAS", "NINF"), ("GA", "NINF"), ("EG", "NINF"),
     ("NPG", "CAS", True), ("NNPG", "NCAS", True), ("GLILD", "NNPG")],
)

# ============================================================================
# EXAMPLE 2 - CD4/CD8 count vs proportion   (fill in after extracting structure)
# EXAMPLE 3 - unusual infection by anatomical site   (fill in likewise)
# ============================================================================


def main():
    render("fig2_granuloma_before", *GRANULOMA_BEFORE)
    render("fig2_granuloma_after", *GRANULOMA_AFTER)
    # render("fig2_tcell_before", *TCELL_BEFORE); ...
    # render("fig2_infection_before", *INFECTION_BEFORE); ...


if __name__ == "__main__":
    main()
