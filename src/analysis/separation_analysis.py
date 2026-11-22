#!/usr/bin/env python3
"""Between-disease separation of individuals in phenotype-similarity space
(Supplementary Figure S4, Supplementary Table S8).

Question: does the curation change how well the three cohorts can be told apart from phenotype
alone, without disease annotations and without a ranking algorithm? Adapted from Campos et al.
2026: minimal term sets, pairwise Lin similarity (symmetric best-match average), multidimensional
scaling on 1 - similarity, PERMANOVA.

Similarity depends on the ontology graph, so each profile is scored in three conditions:

    condition   profile                 graph (ancestors, minimal set, structural IC)
    pre_old     pre-workshop terms      v2024-08-13
    pre_new     pre-workshop terms      v2026-06-23
    post_new    post-workshop terms     v2026-06-23

    comparison  from -> to              what changes
    combined    pre_old -> post_new     terms and hierarchy (the curation as a whole)
    hierarchy   pre_old -> pre_new      graph only (re-parenting, descendant counts)
    terms       pre_new -> post_new     terms only

Settings: observed features only; structural information content
IC(t) = 1 - ln(|desc(t)| + 1) / ln(N) under Phenotypic abnormality (Seco et al. 2004), computed
on the graph of the condition; Lin similarity as primary and Resnik as sensitivity measure.
Separation is summarised by the per-individual silhouette width (cohort as label, the paired
endpoint), the PERMANOVA R2 of the distance matrix and the gap between within- and
between-cohort similarity.

Statistics: similarity between an individual scored on one graph and an individual scored on
the other is not defined, so every statistic is computed within a condition and compared
across conditions:
  * silhouette: paired per individual; Wilcoxon signed-rank per cohort (Holm across cohorts)
    and overall, and a sign-flip permutation of the per-individual differences
  * R2 and gap: leave-one-individual-out jackknife of the difference (SE, 95% CI, z-test)

Outputs (in _work/):
    tables/tableS8a_separation_cohorts.csv         silhouette per cohort and comparison (Lin)
    tables/tableS8b_separation_global.csv          silhouette, R2 and gap per comparison and measure
    intermediate/separation_patient_level.csv      silhouette per individual, condition and measure
    intermediate/separation_term_changes.csv       pre-workshop terms whose ancestors differ between graphs
    figures/figS4_separation_mds.png               multidimensional scaling of the three conditions

Usage, from src/analysis:
    python separation_analysis.py
"""

import math
import pathlib
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import hpotk
from google.protobuf.json_format import Parse
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
from scipy import stats
from scipy.spatial import ConvexHull

HERE       = pathlib.Path(__file__).resolve().parent
FPATH_COH  = HERE.parent / "cohorts"
FPATH_WORK = HERE / "_work"
FPATH_FIGS = FPATH_WORK / "figures"

HPO_OLD_TAG = "v2024-08-13"
HPO_NEW_TAG = "v2026-06-23"
ROOT        = "HP:0000118"           # Phenotypic abnormality
SEED        = 20260813
N_PERM_LABEL = 999                   # PERMANOVA label permutations
N_PERM_FLIP  = 5000                  # sign-flips of the per-individual silhouette difference

# cohort colours, as in the other figures
COHORTS = [("SOCS1", "#003754"), ("APDS1", "#C8A870"), ("NFKB1", "#4E7E96")]
COLOR   = dict(COHORTS)
ORDER   = [c for c, _ in COHORTS]

CONDITIONS = {                       # name -> (profile, graph)
    "pre_old":  ("pre",  "old"),
    "pre_new":  ("pre",  "new"),
    "post_new": ("post", "new"),
}
COMPARISONS = {                      # name -> (from, to, what changes)
    "combined":  ("pre_old", "post_new", "terms and hierarchy"),
    "hierarchy": ("pre_old", "pre_new",  "hierarchy only"),
    "terms":     ("pre_new", "post_new", "terms only"),
}
PRIMARY = "lin"


# --------------------------------------------------------------------------- ontology
def load_hpo(which):
    local = FPATH_WORK / "data" / which / "hp.json"
    if local.exists():
        return hpotk.load_minimal_ontology(str(local))
    tag = HPO_OLD_TAG if which == "old" else HPO_NEW_TAG
    return hpotk.configure_ontology_store().load_minimal_hpo(release=tag)


class Graph:
    """Ancestor sets and structural IC of one HPO release, under Phenotypic abnormality."""

    def __init__(self, hpo):
        self.hpo = hpo
        self.sub_terms = {t.value for t in hpo.graph.get_descendants(ROOT)} | {ROOT}
        self.n_terms = len(self.sub_terms)
        self._anc, self._ndesc = {}, {}

    def ancestors(self, curie):
        """Ancestors incl. self; empty if the term is not under Phenotypic abnormality."""
        a = self._anc.get(curie)
        if a is None:
            if curie in self.sub_terms:
                a = frozenset({t.value for t in self.hpo.graph.get_ancestors(curie)} | {curie})
            else:
                a = frozenset()
            self._anc[curie] = a
        return a

    def ic(self, curie):
        """Seco et al. 2004: 1 - ln(|desc| + 1) / ln(N); None outside the subtree."""
        if curie not in self.sub_terms:
            return None
        v = self._ndesc.get(curie)
        if v is None:
            v = len({t.value for t in self.hpo.graph.get_descendants(curie)}) + 1
            self._ndesc[curie] = v
        return -math.log(v / self.n_terms) / math.log(self.n_terms)

    def minimal_set(self, curies):
        """Drop every term that is an ancestor of another term in the set."""
        s = set(curies)
        redundant = set()
        for c in s:
            redundant |= (self.ancestors(c) - {c}) & s
        return sorted(s - redundant)


def generalize_term(curie, hpo_old, hpo_new):
    # identical to esid4hpo_analysis.ipynb
    if curie in hpo_old:
        return curie
    if curie not in hpo_new:
        return None
    seen, stack = set(), list(hpo_new.graph.get_parents(curie))
    while stack:
        cur = stack.pop()
        if cur.value in seen:
            continue
        seen.add(cur.value)
        if cur in hpo_old:
            return cur.value
        stack.extend(hpo_new.graph.get_parents(cur))
    return None


# --------------------------------------------------------------------------- phenopackets
def observed_curies(pp):
    seen, out = set(), []
    for pf in pp.phenotypic_features:
        if pf.excluded:
            continue
        if pf.type.id not in seen:
            seen.add(pf.type.id)
            out.append(pf.type.id)
    return out


def load_profiles(hpo_old, hpo_new):
    """One row per individual with the post- and the pre-workshop observed term list."""
    rows = []
    for cohort in ORDER:
        ppdir = FPATH_COH / cohort / "phenopackets"
        for fp in sorted(ppdir.glob("*.json")):
            pp = Parse(fp.read_text(), Phenopacket())
            m = re.search(r"PMID[_:]?(\d+)", fp.name)
            post = observed_curies(pp)
            pre = list(dict.fromkeys(g for g in (generalize_term(t, hpo_old, hpo_new) for t in post) if g))
            rows.append(dict(cohort=cohort, phenopacket=fp.name, pmid=m.group(1) if m else "NA",
                             post=post, pre=pre))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- similarity
def similarity_matrix(sets, graph, measure):
    """N x N symmetric best-match average of Lin or Resnik similarity on one graph."""
    terms = sorted({t for s in sets for t in s})
    ic = {t: graph.ic(t) for t in terms}
    idx = {t: i for i, t in enumerate(terms)}
    anc = [graph.ancestors(t) for t in terms]
    ic_anc = {}                                          # IC of every ancestor, cached

    def _ic(a):
        v = ic_anc.get(a)
        if v is None:
            v = graph.ic(a) or 0.0
            ic_anc[a] = v
        return v

    M = np.zeros((len(terms), len(terms)))               # IC of the most informative common ancestor
    for i in range(len(terms)):
        for j in range(i, len(terms)):
            M[i, j] = M[j, i] = max((_ic(a) for a in anc[i] & anc[j]), default=0.0)
    if measure == "lin":
        ic_arr = np.array([ic[t] for t in terms])
        denom = ic_arr[:, None] + ic_arr[None, :]
        T = np.divide(2 * M, denom, out=np.zeros_like(M), where=denom > 0)
    else:                                                # resnik, IC already in 0..1
        T = M
    n = len(sets)
    S = np.eye(n)
    ids = [np.array([idx[t] for t in s], dtype=int) for s in sets]
    for i in range(n):
        for j in range(i + 1, n):
            sub = T[np.ix_(ids[i], ids[j])]
            S[i, j] = S[j, i] = 0.5 * (sub.max(axis=1).mean() + sub.max(axis=0).mean())
    return S


# --------------------------------------------------------------------------- separation metrics
def silhouette(D, labels):
    labels = np.asarray(labels)
    n = len(labels)
    s = np.zeros(n)
    for i in range(n):
        own = (labels == labels[i])
        own[i] = False
        a = D[i, own].mean() if own.any() else 0.0
        b = min(D[i, labels == g].mean() for g in set(labels) if g != labels[i])
        s[i] = 0.0 if max(a, b) == 0 else (b - a) / max(a, b)
    return s


def permanova(D, labels, n_perm=0, rng=None):
    labels = np.asarray(labels)
    n = len(labels)
    groups = sorted(set(labels))
    D2 = D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n

    def ssw(lab):
        tot = 0.0
        for g in groups:
            m = lab == g
            k = m.sum()
            tot += D2[np.ix_(m, m)][np.triu_indices(k, 1)].sum() / k
        return tot

    w = ssw(labels)
    a = sst - w
    F = (a / (len(groups) - 1)) / (w / (n - len(groups)))
    R2 = a / sst
    p = np.nan
    if n_perm:
        cnt = 0
        for _ in range(n_perm):
            w2 = ssw(rng.permutation(labels))
            cnt += ((sst - w2) / (len(groups) - 1)) / (w2 / (n - len(groups))) >= F
        p = (cnt + 1) / (n_perm + 1)
    return F, R2, p


def within_between(S, labels):
    labels = np.asarray(labels)
    iu = np.triu_indices(len(labels), 1)
    same = labels[iu[0]] == labels[iu[1]]
    return S[iu][same].mean(), S[iu][~same].mean()


def global_stats(S, labels):
    D = 1 - S
    _, R2, _ = permanova(D, labels)
    w, b = within_between(S, labels)
    return dict(R2=R2, gap=w - b, within=w, between=b)


def jackknife_delta(S_a, S_b, labels, keys=("R2", "gap")):
    """Leave-one-individual-out jackknife of stat(S_b) - stat(S_a): SE, 95% CI, two-sided z-test p."""
    n = len(labels)
    full = {k: global_stats(S_b, labels)[k] - global_stats(S_a, labels)[k] for k in keys}
    loo = {k: np.zeros(n) for k in keys}
    for i in range(n):
        keep = np.ones(n, bool)
        keep[i] = False
        a = global_stats(S_a[np.ix_(keep, keep)], labels[keep])
        b = global_stats(S_b[np.ix_(keep, keep)], labels[keep])
        for k in keys:
            loo[k][i] = b[k] - a[k]
    out = {}
    for k in keys:
        se = math.sqrt((n - 1) / n * ((loo[k] - loo[k].mean()) ** 2).sum())
        z = full[k] / se if se > 0 else np.nan
        out[k] = dict(delta=full[k], se=se, lo=full[k] - 1.96 * se, hi=full[k] + 1.96 * se,
                      p=2 * stats.norm.sf(abs(z)) if se > 0 else np.nan)
    return out


def holm(ps):
    ps = np.asarray(ps, float)
    order = np.argsort(ps)
    adj = np.empty_like(ps)
    for rank, i in enumerate(order):
        adj[i] = min(1.0, ps[i] * (len(ps) - rank))
    return np.maximum.accumulate(adj[order])[np.argsort(order)]


def cmdscale(D, k=2):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    w, v = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
    order = np.argsort(w)[::-1][:k]
    return v[:, order] * np.sqrt(np.clip(w[order], 0, None))


def procrustes_align(X_ref, X):
    A = X_ref - X_ref.mean(0)
    B = X - X.mean(0)
    U, _, Vt = np.linalg.svd(A.T @ B)
    return B @ (U @ Vt).T


def fmt_p(p):
    if pd.isna(p):
        return ""
    if p >= 1e-3:
        return f"{p:.2g}"
    mant, exp = f"{p:.1e}".split("e-")
    return f"{mant} × 10⁻" + str(int(exp)).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))


# --------------------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(SEED)
    hpo_old, hpo_new = load_hpo("old"), load_hpo("new")
    print("old:", hpo_old.version, "| new:", hpo_new.version)
    graphs = {"old": Graph(hpo_old), "new": Graph(hpo_new)}
    print(f"terms under {ROOT}: old {graphs['old'].n_terms}, new {graphs['new'].n_terms}")

    prof = load_profiles(hpo_old, hpo_new)
    # minimal term set per condition; terms without IC on that graph are dropped
    for cond, (profile, g) in CONDITIONS.items():
        G = graphs[g]
        prof[cond] = prof[profile].map(lambda ts, G=G: G.minimal_set([t for t in ts if G.ic(t) is not None]))
        lost = sorted({t for ts in prof[profile] for t in ts if G.ic(t) is None})
        if lost:
            print(f"[{cond}] {len(lost)} term(s) not under {ROOT} on the {g} graph, dropped: {lost[:5]}")
    # individuals without an observed term cannot be placed in a similarity space
    empty = prof[[c for c in CONDITIONS]].map(len).eq(0).any(axis=1)
    if empty.any():
        print(f"dropped {int(empty.sum())} individual(s) with no observed term:", list(prof.phenopacket[empty]))
        prof = prof[~empty].reset_index(drop=True)
    N = len(prof)
    labels = prof.cohort.values
    print(N, "individuals; mean minimal-set size per condition:")
    print(pd.DataFrame({c: prof[c].map(len).groupby(prof.cohort).mean().round(1) for c in CONDITIONS}).loc[ORDER])

    # how much of the pre-workshop profile is touched by the graph change
    pre_terms = sorted({t for ts in prof.pre for t in ts})
    old_ids = graphs["old"].sub_terms
    moved = {t for t in pre_terms
             if graphs["old"].ancestors(t) & old_ids != graphs["new"].ancestors(t) & old_ids}
    ic_shift = [graphs["new"].ic(t) - graphs["old"].ic(t) for t in pre_terms
                if graphs["old"].ic(t) is not None and graphs["new"].ic(t) is not None]
    print(f"pre-workshop profiles use {len(pre_terms)} distinct terms; {len(moved)} have a different set of "
          f"pre-existing ancestors on the new graph (re-parented, directly or through an ancestor)")
    print(f"  {int(prof.pre.map(lambda ts: any(t in moved for t in ts)).sum())} of {N} individuals carry at least one such term")
    print(f"  structural IC of these terms, new minus old graph: mean {np.mean(ic_shift):+.4f}, "
          f"{int(np.sum(np.abs(ic_shift) > 0.01))} terms shift by more than 0.01")
    (FPATH_WORK / "intermediate").mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(pre_workshop_terms=len(pre_terms), terms_with_changed_ancestors=len(moved),
                       individuals=N,
                       individuals_with_changed_term=int(prof.pre.map(lambda ts: any(t in moved for t in ts)).sum()),
                       mean_ic_shift=round(float(np.mean(ic_shift)), 4))]
                 ).to_csv(FPATH_WORK / "intermediate/separation_term_changes.csv", index=False)

    # similarity per (condition, measure)
    SIM = {}
    for measure in ("lin", "resnik"):
        for cond, (_, g) in CONDITIONS.items():
            SIM[(cond, measure)] = similarity_matrix(list(prof[cond]), graphs[g], measure)

    # per condition: silhouette per individual, global metrics, PERMANOVA label-permutation p
    sil, glob = {}, {}
    for key, S in SIM.items():
        sil[key] = silhouette(1 - S, labels)
        glob[key] = global_stats(S, labels)
        glob[key]["p_permanova"] = permanova(1 - S, labels, N_PERM_LABEL, rng)[2]
        glob[key]["sil_mean"] = sil[key].mean()

    pat_rows = []
    for (cond, measure), s in sil.items():
        for i, r in prof.iterrows():
            pat_rows.append(dict(measure=measure, condition=cond, cohort=r.cohort, phenopacket=r.phenopacket,
                                 pmid=r.pmid, n_terms=len(r[cond]), silhouette=s[i]))
    pd.DataFrame(pat_rows).to_csv(FPATH_WORK / "intermediate/separation_patient_level.csv", index=False)

    # paired comparisons
    coh_rows, glob_rows = [], []
    for measure in ("lin", "resnik"):
        for cmp, (a, b, what) in COMPARISONS.items():
            sa, sb = sil[(a, measure)], sil[(b, measure)]
            d = sb - sa
            p_coh = []
            for g in ORDER:
                m = labels == g
                p_coh.append(stats.wilcoxon(sb[m], sa[m]).pvalue if np.any(d[m] != 0) else np.nan)
            p_holm = holm(p_coh)
            p_all = stats.wilcoxon(sb, sa).pvalue if np.any(d != 0) else np.nan
            signs = rng.choice([-1.0, 1.0], size=(N_PERM_FLIP, N))
            p_flip = (np.sum(np.abs((signs * d).mean(axis=1)) >= abs(d.mean())) + 1) / (N_PERM_FLIP + 1)
            for g in ORDER + ["all"]:
                m = np.ones(N, bool) if g == "all" else labels == g
                coh_rows.append(dict(measure=measure, comparison=cmp, changes=what, cohort=g, n=int(m.sum()),
                                     sil_from=sa[m].mean(), sil_to=sb[m].mean(), sil_delta=d[m].mean(),
                                     improved=int((d[m] > 0).sum()), worse=int((d[m] < 0).sum()),
                                     p_wilcoxon=p_all if g == "all" else p_coh[ORDER.index(g)],
                                     p_holm=np.nan if g == "all" else p_holm[ORDER.index(g)],
                                     p_signflip=p_flip if g == "all" else np.nan))
            jk = jackknife_delta(SIM[(a, measure)], SIM[(b, measure)], labels)
            ga, gb = glob[(a, measure)], glob[(b, measure)]
            glob_rows.append(dict(measure=measure, comparison=cmp, changes=what, **{"from": a, "to": b}, n=N,
                                  sil_from=ga["sil_mean"], sil_to=gb["sil_mean"], sil_delta=d.mean(),
                                  p_sil_wilcoxon=p_all, p_sil_signflip=p_flip,
                                  R2_from=ga["R2"], R2_to=gb["R2"], R2_delta=jk["R2"]["delta"],
                                  R2_delta_lo=jk["R2"]["lo"], R2_delta_hi=jk["R2"]["hi"], p_R2_jackknife=jk["R2"]["p"],
                                  p_permanova_from=ga["p_permanova"], p_permanova_to=gb["p_permanova"],
                                  within_from=ga["within"], within_to=gb["within"],
                                  between_from=ga["between"], between_to=gb["between"],
                                  gap_from=ga["gap"], gap_to=gb["gap"], gap_delta=jk["gap"]["delta"],
                                  gap_delta_lo=jk["gap"]["lo"], gap_delta_hi=jk["gap"]["hi"],
                                  p_gap_jackknife=jk["gap"]["p"]))
    coh = pd.DataFrame(coh_rows)
    gl = pd.DataFrame(glob_rows)

    (FPATH_WORK / "tables").mkdir(parents=True, exist_ok=True)
    t = coh[coh.measure == PRIMARY]
    t = t.assign(**{
        "silhouette (from → to)": t.apply(lambda r: f"{r.sil_from:.3f} → {r.sil_to:.3f}", axis=1),
        "Δ silhouette": t.sil_delta.map(lambda v: f"{v:+.3f}"),
        "improved/worse": t.apply(lambda r: f"{r.improved}/{r.worse}", axis=1),
        "p (Wilcoxon, Holm)": t.apply(lambda r: fmt_p(r.p_wilcoxon if r.cohort == "all" else r.p_holm), axis=1),
    })[["comparison", "changes", "cohort", "n", "silhouette (from → to)", "Δ silhouette", "improved/worse",
        "p (Wilcoxon, Holm)"]]
    t.to_csv(FPATH_WORK / "tables/tableS8a_separation_cohorts.csv", index=False)
    gl.round(4).to_csv(FPATH_WORK / "tables/tableS8b_separation_global.csv", index=False)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print("\nTable S8a (Lin, per cohort)\n", t.to_string(index=False))
    print("\nTable S8b (global)\n", gl[["measure", "comparison", "sil_from", "sil_to", "p_sil_wilcoxon", "p_sil_signflip",
                                         "R2_from", "R2_to", "R2_delta", "R2_delta_lo", "R2_delta_hi", "p_R2_jackknife",
                                         "gap_from", "gap_to", "gap_delta", "p_gap_jackknife"]].round(4).to_string(index=False))

    # ---- figure: MDS of the three conditions (primary measure), aligned to the first panel
    C_INK, C_MUTED = "#003754", "#5f7078"
    mpl.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": None,
        "figure.facecolor": "white", "axes.facecolor": "white", "font.size": 10.5,
        "font.family": "DejaVu Sans", "text.color": C_INK, "axes.titlesize": 10.5,
        "axes.labelsize": 10, "axes.labelcolor": C_MUTED, "axes.edgecolor": "#c4c9ce",
        "axes.spines.top": False, "axes.spines.right": False, "xtick.color": C_MUTED,
        "ytick.color": C_INK, "legend.frameon": False, "legend.fontsize": 10})
    titles = {"pre_old":  f"pre-workshop terms\n{HPO_OLD_TAG} graph",
              "pre_new":  f"pre-workshop terms\n{HPO_NEW_TAG} graph",
              "post_new": f"post-workshop terms\n{HPO_NEW_TAG} graph"}
    X_ref = cmdscale(1 - SIM[("pre_old", PRIMARY)])
    X_ref = X_ref - X_ref.mean(0)
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.4), sharex=True, sharey=True)
    for ax, letter, cond in zip(axes, ("a)", "b)", "c)"), CONDITIONS):
        X = X_ref if cond == "pre_old" else procrustes_align(X_ref, cmdscale(1 - SIM[(cond, PRIMARY)]))
        for g in ORDER:
            m = labels == g
            P = X[m]
            ax.scatter(P[:, 0], P[:, 1], s=20, color=COLOR[g], alpha=0.75, lw=0, label=f"{g} (n = {m.sum()})")
            if m.sum() >= 3:
                poly = P[ConvexHull(P).vertices]
                ax.fill(poly[:, 0], poly[:, 1], color=COLOR[g], alpha=0.10, lw=0)
                ax.plot(np.r_[poly[:, 0], poly[0, 0]], np.r_[poly[:, 1], poly[0, 1]], color=COLOR[g], lw=0.8, alpha=0.6)
            ax.scatter(*P.mean(0), marker="X", s=120, color=COLOR[g], edgecolor="white", lw=1.2, zorder=5)
        st = glob[(cond, PRIMARY)]
        ax.set_title(f"{titles[cond]}\nR² = {st['R2']:.3f}, silhouette = {st['sil_mean']:.3f}")
        ax.set_xlabel("MDS dimension 1")
        ax.set_aspect("equal", adjustable="box")
        ax.text(-0.02, 1.2, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", color=C_INK,
                ha="right", va="bottom")
    axes[0].set_ylabel("MDS dimension 2")
    handles, labs = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labs, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    FPATH_FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FPATH_FIGS / "figS4_separation_mds.png")
    print("saved figS4_separation_mds.png")


if __name__ == "__main__":
    main()
