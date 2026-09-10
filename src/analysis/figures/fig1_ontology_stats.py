#!/usr/bin/env python3
"""
Figure 1 - the ESID4HPO curation effort and its output.

Every term under "Abnormality of the immune system" (HP:0002715) is treated as
in scope. Per subbranch:

  panel A  Terms added        existing vs newly created terms
  panel B  Existing improved  re-parented / redefined / obsoleted
  panel C  Synonyms added     alternative labels for findability

is_a edge churn and mean depth are written to the summary CSV (reported in text
/ supplement rather than plotted).

Inputs are the two release hp.json files of the LIRICAL data bundles
(src/analysis/_work/data/{old,new}/hp.json); outputs go to
src/analysis/_work/figures/.

USAGE (from src/analysis/figures)
  python fig1_ontology_stats.py --list-anchors
  python fig1_ontology_stats.py
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import hpotk

# =============================================================================
# CONFIG
# =============================================================================

ANALYSIS = pathlib.Path(__file__).resolve().parents[1]    # src/analysis
HP_OLD = ANALYSIS / "_work" / "data" / "old" / "hp.json"     # v2024-08-13
HP_NEW = ANALYSIS / "_work" / "data" / "new" / "hp.json"     # v2026-06-23
OUTDIR = ANALYSIS / "_work" / "figures"

OLD_TAG = "v2024-08-13"
NEW_TAG = "v2026-06-23"

ROOTS = ["HP:0002715"]        # Abnormality of the immune system

GROUPS = {
    "Flow cytometry /\nimmunophenotyping": ["HP:0010987"],
    "Antibodies /\nhumoral":               ["HP:0005368"],
    "Infections":                          ["HP:0032101"],
}
OTHER_LABEL = "Other immune"

# ---- palette: BIH brand colours, lightened -----------------------------------
C_EXIST = "#dde4e8"   # pale tint of BIH navy    (existing terms)
C_NEW   = "#4e7e96"   # BIH navy/teal, softened  (new terms)
C_REPAR = "#c8a870"   # BIH gold, softened       (re-parented)
C_DEF   = "#6cbcc5"   # BIH teal, softened       (redefined)
C_OBS   = "#c9cdd1"   # neutral grey             (obsoleted)
C_SYN   = "#a6a3cd"   # BIH purple, softened     (synonyms)
INK, MUTED = "#003754", "#5f7078"   # BIH navy for text; muted blue-grey

FIGSIZE = (13.6, 4.6)
BARH = 0.58


# =============================================================================
# loading
# =============================================================================

def load_graphs():
    for p in (HP_OLD, HP_NEW):
        if not p.exists():
            sys.exit(f"ERROR: {p} not found (expected the LIRICAL data bundles).")
    old = hpotk.load_minimal_ontology(str(HP_OLD))
    new = hpotk.load_minimal_ontology(str(HP_NEW))
    print(f"loaded graphs: old={old.version}  new={new.version}")
    return old, new


def _curie(iri):
    return iri.rsplit("/", 1)[-1].replace("_", ":")


def load_metadata(path):
    with path.open() as fh:
        doc = json.load(fh)
    out = {}
    for graph in doc.get("graphs", []):
        for node in graph.get("nodes", []):
            nid = node.get("id", "")
            if "/HP_" not in nid:
                continue
            meta = node.get("meta", {}) or {}
            out[_curie(nid)] = dict(
                label=node.get("lbl"),
                definition=(meta.get("definition") or {}).get("val"),
                synonyms={s.get("val") for s in meta.get("synonyms", []) if s.get("val")},
                deprecated=bool(meta.get("deprecated", False)),
            )
    return out


# =============================================================================
# anchors
# =============================================================================

def list_anchors(new):
    for root in ROOTS:
        if root not in new:
            print(f"\n!! {root} not present in {NEW_TAG}")
            continue
        print(f"\n=== children of {root} ({new.get_term_name(root)}) ===")
        rows = [(sum(1 for _ in new.graph.get_descendants(ch, include_source=True)),
                 ch.value, new.get_term_name(ch))
                for ch in new.graph.get_children(root)]
        for n, curie, label in sorted(rows, reverse=True):
            print(f"  {curie}  {n:5d} descendants   {label}")


# =============================================================================
# term table
# =============================================================================

def descendants(o, curie):
    return {t.value for t in o.graph.get_descendants(curie, include_source=True)} \
        if curie in o else set()


def parents(o, curie):
    return frozenset(p.value for p in o.graph.get_parents(curie))


def depth_map(o, roots):
    dist, queue = {}, collections.deque()
    for r in roots:
        if r in o:
            dist[r] = 0
            queue.append(r)
    while queue:
        cur = queue.popleft()
        for ch in o.graph.get_children(cur):
            if ch.value not in dist:
                dist[ch.value] = dist[cur] + 1
                queue.append(ch.value)
    return dist


def build_term_table(old, new, meta_old, meta_new):
    domain_new = set().union(*(descendants(new, r) for r in ROOTS))
    domain_old = set().union(*(descendants(old, r) for r in ROOTS))
    members = {g: (set().union(*(descendants(new, a) for a in anchors)) if anchors else set())
               for g, anchors in GROUPS.items()}

    def classify(curie):
        for g in GROUPS:
            if curie in members[g]:
                return g
        return OTHER_LABEL

    d_old, d_new = depth_map(old, ROOTS), depth_map(new, ROOTS)
    rows = []
    for curie in sorted(domain_new):
        m_new, m_old = meta_new.get(curie, {}), meta_old.get(curie, {})
        is_new = curie not in old
        if is_new:
            rel_added = rel_removed = 0
        else:
            p_old, p_new = parents(old, curie), parents(new, curie)
            rel_added, rel_removed = len(p_new - p_old), len(p_old - p_new)
        rows.append(dict(
            curie=curie, label=m_new.get("label") or new.get_term_name(curie),
            group=classify(curie), in_old=not is_new, new_term=is_new,
            reparented=(not is_new) and (rel_added or rel_removed) > 0,
            rel_added=rel_added, rel_removed=rel_removed,
            def_changed=(not is_new) and m_old.get("definition") != m_new.get("definition"),
            syn_added=len((m_new.get("synonyms") or set()) - (m_old.get("synonyms") or set())),
            depth_new=d_new.get(curie, np.nan), depth_old=d_old.get(curie, np.nan),
            obsoleted=False))
    for curie in sorted(domain_old - domain_new):
        rows.append(dict(
            curie=curie, label=meta_old.get(curie, {}).get("label"),
            group=OTHER_LABEL, in_old=True, new_term=False, reparented=False,
            rel_added=0, rel_removed=0, def_changed=False, syn_added=0,
            depth_new=np.nan, depth_old=d_old.get(curie, np.nan), obsoleted=True))
    return pd.DataFrame(rows)


def summarise(df):
    order = list(GROUPS) + [OTHER_LABEL]
    live = df[~df.obsoleted]
    s = (live.groupby("group")
             .agg(terms_after=("curie", "size"), new_terms=("new_term", "sum"),
                  reparented=("reparented", "sum"), rel_added=("rel_added", "sum"),
                  rel_removed=("rel_removed", "sum"), def_changed=("def_changed", "sum"),
                  syn_added=("syn_added", "sum"), mean_depth_after=("depth_new", "mean"))
             .reindex(order))
    s["terms_before"] = live[live.in_old].groupby("group").size().reindex(order)
    s["obsoleted"] = df[df.obsoleted].groupby("group").size().reindex(order)
    s = s.fillna(0)
    ints = [c for c in s.columns if not c.startswith("mean_depth")]
    s[ints] = s[ints].astype(int)
    s["pct_growth"] = np.where(s.terms_before > 0,
                               (s.terms_after - s.terms_before) / s.terms_before * 100,
                               np.nan).round(1)
    return s


# =============================================================================
# figure
# =============================================================================

def _style():
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.size": 10.5, "font.family": "DejaVu Sans", "text.color": INK,
        "axes.titlesize": 11.5, "axes.titlepad": 8,
        "axes.labelsize": 10, "axes.labelcolor": MUTED,
        "axes.edgecolor": "#c4c9ce", "axes.linewidth": 1.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": INK,
        "xtick.labelsize": 9.5, "ytick.labelsize": 10,
        "legend.frameon": False, "legend.fontsize": 9,
    })


def _stack(ax, y, segments):
    """segments = [(values, color, label), ...]; returns per-row totals."""
    left = np.zeros(len(y))
    for vals, color, label in segments:
        vals = np.asarray(vals, float)
        ax.barh(y, vals, BARH, left=left, color=color,
                edgecolor="white", linewidth=0.9, label=label, zorder=3)
        left += vals
    return left


def _clean(ax, xmax, n):
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.7, n - 0.3)                      # margin so bottom bar clears the axis
    ax.grid(axis="x", color="#000", alpha=.06, lw=.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    ax.margins(y=0)


def _labels(ax, y, xpos, values, xmax, color, weight="normal", plus=False):
    for yi, xp, v in zip(y, xpos, values):
        if v:
            ax.text(xp + xmax * .015, yi, (f"+{int(v)}" if plus else f"{int(v)}"),
                    va="center", fontsize=9, color=color, fontweight=weight)


def _legend_below(ax, colors, labels, ncol):
    handles = [mpl.patches.Patch(facecolor=c, edgecolor="none") for c in colors]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=ncol, handlelength=1.0, handleheight=1.0, columnspacing=1.3,
              borderaxespad=0)


def make_figure(summary, outdir):
    _style()
    order = list(summary.index)[::-1]
    s = summary.loc[order]
    y = np.arange(len(s))
    n = len(s)

    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=FIGSIZE, sharey=True, gridspec_kw=dict(wspace=0.14))

    # ---- A: terms added ----------------------------------------------------
    top = _stack(axA, y, [(s.terms_before.values, C_EXIST, "existing"),
                          (s.new_terms.values, C_NEW, "new")])
    xmaxA = max(top.max(), 1) * 1.16
    _clean(axA, xmaxA, n)
    _labels(axA, y, top, s.new_terms.values, xmaxA, C_NEW, weight="bold", plus=True)
    axA.set_yticks(y)
    axA.set_yticklabels(s.index)
    axA.set_xlabel("Terms in branch")
    axA.set_title("I) Terms added", loc="left", fontweight="bold", color=INK)
    _legend_below(axA, [C_EXIST, C_NEW], ["existing", "new"], ncol=2)

    # ---- B: existing terms improved ---------------------------------------
    top = _stack(axB, y, [(s.reparented.values, C_REPAR, "re-parented"),
                          (s.def_changed.values, C_DEF, "redefined"),
                          (s.obsoleted.values, C_OBS, "obsoleted")])
    xmaxB = max(top.max(), 1) * 1.18
    _clean(axB, xmaxB, n)
    _labels(axB, y, top, top, xmaxB, MUTED)
    axB.set_xlabel("Existing terms affected")
    axB.set_title("II) Existing terms improved", loc="left", fontweight="bold", color=INK)
    _legend_below(axB, [C_REPAR, C_DEF, C_OBS],
                  ["re-parented", "redefined", "obsoleted"], ncol=3)

    # ---- C: synonyms added -------------------------------------------------
    _stack(axC, y, [(s.syn_added.values, C_SYN, "synonyms")])
    xmaxC = max(s.syn_added.max(), 1) * 1.16
    _clean(axC, xmaxC, n)
    _labels(axC, y, s.syn_added.values, s.syn_added.values, xmaxC, MUTED)
    axC.set_xlabel("Synonyms added")
    axC.set_title("III) Synonyms added", loc="left", fontweight="bold", color=INK)
    _legend_below(axC, [C_SYN], ["synonyms"], ncol=1)

    fig.subplots_adjust(left=0.17, right=0.985, top=0.87, bottom=0.22)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "figure1_panels.png")
    fig.savefig(outdir / "figure1_panels.pdf")
    print(f"figure -> {outdir / 'figure1_panels.png'}")


# =============================================================================
# main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list-anchors", action="store_true")
    args = ap.parse_args()

    old, new = load_graphs()
    if args.list_anchors:
        list_anchors(new)
        return

    meta_new, meta_old = load_metadata(HP_NEW), load_metadata(HP_OLD)
    df = build_term_table(old, new, meta_old, meta_new)
    summary = summarise(df)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTDIR / "fig1_term_level.csv", index=False)
    summary.to_csv(OUTDIR / "fig1_summary.csv")

    with pd.option_context("display.width", 220, "display.max_columns", 50):
        print("\n", summary, "\n")
    print(f"total terms in domain: {(~df.obsoleted).sum()} "
          f"(new since {OLD_TAG}: {int(df.new_term.sum())}); "
          f"is_a edges added={int(df.rel_added.sum())}, removed={int(df.rel_removed.sum())}, "
          f"synonyms added={int(df.syn_added.sum())}")

    make_figure(summary, OUTDIR)


if __name__ == "__main__":
    main()
