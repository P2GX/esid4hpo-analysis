#!/usr/bin/env python3
"""
Figure 2 - the kinds of curation applied to the immune branch, one worked example each.

    a) anatomy of an HPO term, and how an existing term was revised
    b) adding a new term to the existing hierarchy
    c) merging / removing a term
    d) restructuring a branch (memory CD4+ T cells)

All labels, definitions, synonyms and is_a relationships are read from the two releases
(_work/data/old/hp.json = v2024-08-13, _work/data/new/hp.json = v2026-06-23); the script
stops if a term or relationship it draws does not exist in the release it is drawn for.
The logical definition in panel a is not part of hp.json; it is taken from hp-base.obo
(identical in both releases) and written out below.

Visual grammar
    blue fill       new term
    teal border     relabelled term
    grey            existing term, unchanged
    dashed          obsoleted term
    white           parent shown for orientation
    ochre arrow     added parent relationship
    grey arrow      existing is_a relationship (child -> parent)

Run:    python fig2_curation_examples.py        (from src/analysis/figures)
Output: src/analysis/_work/figures/fig2_curation_examples.png
"""

import json
import pathlib
import textwrap

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = pathlib.Path(__file__).resolve().parent
WORK = HERE.parent / "_work"
OUT = WORK / "figures"
OLD_TAG, NEW_TAG = "v2024-08-13", "v2026-06-23"

NAVY, BLUE, GREY = "#003754", "#4E7E96", "#E7EAEC"
AMBER, TEAL = "#C8A870", "#6CBCC5"
INK, MUT, BORD, ARROW = "#2B2F33", "#5F7078", "#C4C9CE", "#9AA3AD"

# text sizes (pt); the canvas is 10 in wide
FS = dict(box=8.2, id=7.2, body=8.6, head=10.5, title=11.5, letter=14, note=8.2, legend=9)
CHAR_IN = 0.0078          # approximate width of one character in inches per pt of font size


# --------------------------------------------------------------------------- ontology
def load(which):
    g = json.load(open(WORK / "data" / which / "hp.json"))["graphs"][0]
    cur = lambda iri: iri.rsplit("/", 1)[-1].replace("_", ":")
    nodes = {cur(n["id"]): n for n in g["nodes"] if "HP_" in n["id"]}
    parents = {}
    for e in g["edges"]:
        if e["pred"] == "is_a" and "HP_" in e["sub"] and "HP_" in e["obj"]:
            parents.setdefault(cur(e["sub"]), set()).add(cur(e["obj"]))
    return nodes, parents


ONT = {"old": load("old"), "new": load("new")}


def label(rel, t):
    return ONT[rel][0][t]["lbl"]


def definition(rel, t):
    return ONT[rel][0][t]["meta"]["definition"]["val"]


def synonyms(rel, t):
    return sorted((s["val"] for s in ONT[rel][0][t].get("meta", {}).get("synonyms", [])), key=str.lower)


def is_obsolete(rel, t):
    return bool(ONT[rel][0][t].get("meta", {}).get("deprecated"))


def replaced_by(rel, t):
    for b in ONT[rel][0][t].get("meta", {}).get("basicPropertyValues", []):
        if b["pred"].endswith("IAO_0100001"):
            return b["val"]
    return None


def check_edge(rel, child, parent):
    assert parent in ONT[rel][1].get(child, ()), f"{rel}: {child} is_a {parent} not in the release"


def kind_of(t):
    """How a term of the new release relates to the old one."""
    if t not in ONT["old"][0]:
        return "new"
    return "relab" if label("old", t) != label("new", t) else "exist"


# --------------------------------------------------------------------------- drawing
fig = plt.figure(figsize=(10, 13.3))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 10)
ax.set_ylim(13.3, 0)            # inches, y pointing down
ax.axis("off")
mpl.rcParams["font.family"] = "DejaVu Sans"

STYLE = {  # fill, edge, linewidth, text colour, bold, linestyle
    "ctx":   ("white", BORD, 1.0, MUT, False, "-"),          # parent shown for orientation
    "exist": (GREY, BORD, 1.0, INK, False, "-"),
    "new":   (BLUE, BLUE, 1.0, "white", True, "-"),
    "relab": (GREY, TEAL, 2.2, INK, False, "-"),
    "obs":   ("white", ARROW, 1.2, MUT, False, (0, (4, 3))),
}


def wrap(text, width_in, fs):
    return "\n".join(textwrap.wrap(text, max(8, int(width_in / (CHAR_IN * fs))), break_long_words=False))


def box(x, y, w, h, term, rel, kind=None, text=None):
    """Term box with wrapped label and the ID below it; returns (centre x, top y, bottom y)."""
    kind = kind or ("exist" if rel == "old" else kind_of(term))
    fill, edge, lw, tc, bold, ls = STYLE[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                fc=fill, ec=edge, lw=lw, ls=ls, zorder=2))
    lab = text or label(rel, term)
    ax.text(x + w / 2, y + h / 2 - 0.075, wrap(lab, (w - 0.12) * (0.94 if bold else 1), FS["box"]), ha="center", va="center",
            fontsize=FS["box"], color=tc, fontweight="bold" if bold else "normal", linespacing=1.12, zorder=3)
    ax.text(x + w / 2, y + h - 0.045, term, ha="center", va="bottom", fontsize=FS["id"],
            color=tc if kind == "new" else MUT, zorder=3)
    return x + w / 2, y, y + h


def arrow(child, parent, added=False, dashed=False):
    (cx, ct, _), (px, _, pb) = child, parent
    ax.add_patch(FancyArrowPatch((cx, ct), (px, pb), arrowstyle="-|>", mutation_scale=11,
                                 color=AMBER if added else ARROW, lw=1.9 if added else 1.1,
                                 ls=(0, (4, 3)) if dashed else "-", shrinkA=0, shrinkB=1, zorder=1))


def edge(rel, c_id, p_id, c_pos, p_pos):
    """is_a arrow, verified against the release; ochre if the relationship is new."""
    check_edge(rel, c_id, p_id)
    added = rel == "new" and (c_id not in ONT["old"][0] or p_id not in ONT["old"][1].get(c_id, ()))
    arrow(c_pos, p_pos, added=added)


def panel_title(letter, x, y, title):
    ax.text(x, y, letter, fontsize=FS["letter"], fontweight="bold", color=NAVY, va="center")
    ax.text(x + 0.38, y, title, fontsize=FS["title"], fontweight="bold", color=NAVY, va="center")


def head(x, y, w, text):
    ax.text(x + w / 2, y, text, fontsize=FS["head"], fontweight="bold", color=MUT, ha="center", va="center")


def note(x, y, w, text, italic=True, ha="left"):
    ax.text(x if ha == "left" else x + w / 2, y, wrap(text, w, FS["note"]), fontsize=FS["note"], color=MUT,
            style="italic" if italic else "normal", va="top", ha=ha, linespacing=1.25)


# =========================================================================== a) anatomy of a term
A_TERM = "HP:0010976"
# from hp-base.obo, identical in v2024-08-13 and v2026-06-23:
#   intersection_of: PATO:0001997 / RO:0000052 CL:0000236 / RO:0002573 PATO:0000460
A_LOGICAL = ("'decreased amount' (PATO:0001997) and 'characteristic of' some 'B cell' (CL:0000236) "
             "and 'has modifier' some 'abnormal' (PATO:0000460)")
panel_title("a)", 0.3, 0.3, "Anatomy of an HPO term, and how an existing term was revised")

CARD_W, LAB_W, CARD_Y = 4.55, 0.95, 0.84
FS_A, LINE = 8.2, 0.15           # body size of the two term cards, line pitch of the synonym list
syn_old, syn_new = synonyms("old", A_TERM), synonyms("new", A_TERM)
assert label("old", A_TERM) != label("new", A_TERM) and definition("old", A_TERM) != definition("new", A_TERM)
ROWS = [("Label", 0.42), ("ID", 0.36), ("Definition", 0.66), ("Synonyms", LINE * len(syn_new) + 0.14),
        ("Logical\ndefinition", 0.62)]
# what happened to each field of the term, shown under the field name in the "after" card
CHANGE = {"Label": ("relabelled", TEAL), "ID": ("unchanged", MUT), "Definition": ("redefined", TEAL),
          "Synonyms": (f"{len(set(syn_new) - set(syn_old))} added", BLUE), "Logical\ndefinition": ("unchanged", MUT)}
CARD_H = sum(h for _, h in ROWS) + 0.06

for rel, cx, tag in (("old", 0.3, f"Before   ({OLD_TAG})"), ("new", 5.15, f"After   ({NEW_TAG})")):
    head(cx, 0.66, CARD_W, tag)
    after = rel == "new"
    ax.add_patch(FancyBboxPatch((cx, CARD_Y), CARD_W, CARD_H, boxstyle="round,pad=0,rounding_size=0.06",
                                fc="white", ec=BORD, lw=1.0, zorder=1))
    y = CARD_Y + 0.03
    for i, (row, h) in enumerate(ROWS):
        if i:
            ax.plot([cx + 0.1, cx + CARD_W - 0.1], [y, y], color=BORD, lw=0.7, zorder=2)
        n_lines = row.count("\n") + 1
        ax.text(cx + 0.18, y + 0.09, row, fontsize=FS_A, fontweight="bold", color=MUT, va="top", linespacing=1.15)
        if after:
            txt, col = CHANGE[row]
            ax.text(cx + 0.18, y + 0.1 + 0.135 * n_lines, txt, fontsize=FS["id"], color=col, va="top",
                    fontweight="normal" if col == MUT else "bold")
        tx, tw = cx + 0.18 + LAB_W, CARD_W - LAB_W - 0.32
        if row == "Label":
            ax.add_patch(FancyBboxPatch((tx - 0.08, y + 0.07), tw + 0.16, h - 0.14,
                                        boxstyle="round,pad=0,rounding_size=0.05", fc=GREY,
                                        ec=TEAL if after else BORD, lw=2.0 if after else 1.0, zorder=2))
            ax.text(tx, y + h / 2, label(rel, A_TERM), fontsize=9.6, fontweight="bold", color=INK, va="center",
                    zorder=3)
        elif row == "ID":
            ax.text(tx, y + 0.09, A_TERM, fontsize=FS_A, color=INK, va="top")
        elif row == "Definition":
            ax.text(tx, y + 0.09, wrap(definition(rel, A_TERM), tw, FS_A), fontsize=FS_A, color=INK, va="top",
                    linespacing=1.2)
        elif row == "Synonyms":
            for k, syn in enumerate(syn_new if after else syn_old):
                added = after and syn not in syn_old
                ax.text(tx, y + 0.09 + LINE * k, ("+ " if added else "") + syn, fontsize=FS_A, va="top",
                        color=BLUE if added else INK, fontweight="bold" if added else "normal")
        else:
            ax.text(tx, y + 0.09, wrap(A_LOGICAL, tw, FS_A), fontsize=FS_A, color=INK, va="top", linespacing=1.2)
        y += h
A_END = CARD_Y + CARD_H

# =========================================================================== b) adding a new term
B_Y = A_END + 0.45
B_X, B_W = 0.3, 4.0
B_NEW, B_PARENT, B_TOP = "HP:5210281", "HP:0020085", "HP:0032101"
B_SIBS = ["HP:0020088", "HP:0020089"]
panel_title("b)", B_X, B_Y, "Adding a new term")
head(B_X, B_Y + 0.4, B_W, f"After   ({NEW_TAG})")
assert B_NEW not in ONT["old"][0] and B_NEW in ONT["new"][0]
bw, bh, gap = 1.24, 0.62, 0.14
y0 = B_Y + 0.65
p_top = box(B_X + (B_W - 2.6) / 2, y0, 2.6, 0.5, B_TOP, "new", kind="ctx")
p_par = box(B_X + (B_W - 2.6) / 2, y0 + 0.82, 2.6, 0.5, B_PARENT, "new")
edge("new", B_PARENT, B_TOP, p_par, p_top)
for i, t in enumerate(B_SIBS + [B_NEW]):
    pos = box(B_X + i * (bw + gap), y0 + 1.66, bw, bh, t, "new")
    edge("new", t, B_PARENT, pos, p_par)
n_children = sum(1 for c, ps in ONT["new"][1].items() if B_PARENT in ps)
note(B_X, y0 + 2.42, B_W, f"{n_children - 3} further child terms not shown.")
note(B_X, y0 + 2.68, B_W, "Definition: " + definition("new", B_NEW))
note(B_X, y0 + 3.3, B_W, "Synonyms: " + "; ".join(synonyms("new", B_NEW)))

# =========================================================================== c) merging / removing
C_X, C_W = 4.6, 5.1
C_PARENT, C_OBS, C_TOP = "HP:0002720", "HP:0011837", "HP:0410240"
C_NEW = ["HP:0003460", "HP:0033024", "HP:5210005"]
panel_title("c)", C_X, B_Y, "Merging or removing a term")
assert is_obsolete("new", C_OBS) and replaced_by("new", C_OBS) == C_PARENT
cw, ch_, cg = 1.2, 0.72, 0.1
# only the release after curation is drawn; the obsoleted term is kept as a dashed box
check_edge("old", C_OBS, C_PARENT)
ya = B_Y + 0.65 + 0.82
head(C_X, B_Y + 0.4, C_W, f"After   ({NEW_TAG})")
c_top = box(C_X + (C_W - 2.8) / 2, ya - 0.82, 2.8, 0.5, C_TOP, "new", kind="ctx")
cp2 = box(C_X + (C_W - 2.8) / 2, ya, 2.8, 0.5, C_PARENT, "new")
edge("new", C_PARENT, C_TOP, cp2, c_top)
ghost = box(C_X, ya + 0.8, cw, ch_, C_OBS, "new", kind="obs", text=label("old", C_OBS))
arrow(ghost, cp2, dashed=True)
ax.text(C_X + cw / 2 - 0.02, ya + 0.63, "merged", fontsize=FS["id"], color=MUT, style="italic", ha="right")
for i, t in enumerate(C_NEW):
    pos = box(C_X + (i + 1) * (cw + cg), ya + 0.8, cw, ch_, t, "new")
    edge("new", t, C_PARENT, pos, cp2)
note(C_X, ya + 1.74, C_W,
     f"{label('old', C_OBS)} ({C_OBS}) had an ambiguous label. It was obsoleted and merged into its parent, which "
     f"keeps {C_OBS} as an alternative ID, so existing annotations remain valid. The child terms now follow the "
     "chronic / transient pattern.")
BC_END = y0 + 3.7

# =========================================================================== d) restructuring a branch
D_Y = BC_END + 0.35
panel_title("d)", 0.3, D_Y, "Restructuring a branch: the memory CD4+ T cell hierarchy")
# before
DB_X, DB_W = 0.3, 3.9
head(DB_X, D_Y + 0.42, DB_W, f"Before   ({OLD_TAG})")
mh = 1.12                      # the pre-curation labels are long and need five lines
yt, ym = D_Y + 0.7, D_Y + 1.68
yl = ym + mh + 0.38
ax.plot([4.4, 4.4], [D_Y + 0.3, yl + mh], color=BORD, lw=0.9)
mw = (DB_W - 0.24) / 3
top = box(DB_X + 0.35, yt, DB_W - 0.7, 0.6, "HP:0410379", "old")
mids = {}
for i, t in enumerate(["HP:0025620", "HP:0410381", "HP:0410391"]):
    mids[t] = box(DB_X + i * (mw + 0.12), ym, mw, mh, t, "old")
    edge("old", t, "HP:0410379", mids[t], top)
lw_ = 1.55
for i, t in enumerate(["HP:0410393", "HP:0410394"]):
    pos = box(DB_X + DB_W - 2 * lw_ - 0.12 + i * (lw_ + 0.12), yl, lw_, mh, t, "old")
    edge("old", t, "HP:0410391", pos, mids["HP:0410391"])
# after
DA_X, DA_W = 4.6, 5.1
head(DA_X, D_Y + 0.42, DA_W, f"After   ({NEW_TAG})")
ag = 0.1
aw = (DA_W - 3 * ag) / 4
top = box(DA_X + 1.05, yt, DA_W - 2.1, 0.6, "HP:0410379", "new")
mid_ids = ["HP:5210410", "HP:0410386", "HP:0410391", "HP:5210409"]
mids = {}
for i, t in enumerate(mid_ids):
    mids[t] = box(DA_X + i * (aw + ag), ym, aw, mh, t, "new")
    edge("new", t, "HP:0410379", mids[t], top)
low = {"HP:0410388": ("HP:5210410", "HP:0410386"), "HP:0410393": ("HP:5210410", "HP:0410391"),
       "HP:0025624": ("HP:0410386", "HP:5210409"), "HP:0025625": ("HP:0410391", "HP:5210409")}
for i, (t, ps) in enumerate(low.items()):
    pos = box(DA_X + i * (aw + ag), yl, aw, mh, t, "new")
    for p in ps:
        edge("new", t, p, pos, mids[p])
D_END = yl + mh

# =========================================================================== legend
ly = D_END + 0.42
lx = 0.3
for kind, txt in (("new", "new term"), ("relab", "relabelled term"), ("exist", "existing term"), ("obs", "obsoleted term")):
    fill, edge_c, lw, _, _, ls = STYLE[kind]
    ax.add_patch(FancyBboxPatch((lx, ly - 0.09), 0.3, 0.18, boxstyle="round,pad=0,rounding_size=0.03",
                                fc=fill, ec=edge_c, lw=min(lw, 1.8), ls=ls))
    ax.text(lx + 0.38, ly, txt, fontsize=FS["legend"], color=INK, va="center")
    lx += 0.38 + 0.072 * len(txt) + 0.2
for col, lw, txt in ((AMBER, 1.9, "added parent"), (ARROW, 1.1, "existing parent")):
    ax.add_patch(FancyArrowPatch((lx, ly), (lx + 0.36, ly), arrowstyle="-|>", mutation_scale=11, color=col, lw=lw))
    ax.text(lx + 0.44, ly, txt, fontsize=FS["legend"], color=INK, va="center")
    lx += 0.44 + 0.072 * len(txt) + 0.2
ax.text(5.0, ly + 0.3, "Arrows point from child to parent (is_a). Terms, labels and relationships are read from "
        f"the HPO releases {OLD_TAG} and {NEW_TAG}.", fontsize=FS["note"], color=MUT, ha="center", va="center")

fig.set_size_inches(10, ly + 0.6)
ax.set_ylim(ly + 0.6, 0)
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "fig2_curation_examples.png", dpi=300, facecolor="white")
print("saved", OUT / "fig2_curation_examples.png", f"(10 x {ly + 0.6:.1f} in)")
