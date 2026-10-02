#!/usr/bin/env python3
"""Figure 5 and Supplementary Figure: ESID4HPO terms recovered from routine clinical text (Southampton).

Figure 5 (fig5_ehr.*), one row per data source in both panels:
  a) per-patient mean IC per term, pre- vs post-workshop release (paired boxes),
     with the mean change and the Wilcoxon p above each pair;
  b) mean change per patient (post - pre) with 95 % bootstrap CI for the three summary
     measures: immune-branch terms per patient, mean IC per term, summed IC.
Supplementary Figure (figS_ehr_terms.*):
  per data source, the terms recoverable only with the post-workshop release, ranked by
  the number of patients, navy = term created by ESID4HPO, gold = pre-existing term that
  became retrievable through re-parenting into the immune branch or relabelling.

Inputs: the summary table (default: Alex's numbers, _work/ehr_ic_summary_sa.csv) for N,
means and p, and the per-patient / per-term tables written by ehr_ic.py for the boxes,
the bootstrap CIs and the term lists. Style and palette follow Figures 1, 3 and 4 (BIH
navy = post-workshop, blue-grey = pre-workshop); legends below the panels; line icons
mark the data source of each row (--no-icons to drop them). A laboratory-value unit can
be added later as one more entry in UNIT_ORDER / UNIT_LABEL and an icon kind.

Usage, from src/analysis:
    python figures/fig5_ehr.py                 # Figure 5 and the supplementary figure
    python figures/fig5_ehr.py --only main     # or: --only supp
Outputs: _work/figures/figure5/fig5_ehr.{png,pdf,svg}, figS_ehr_terms.{png,pdf,svg}
"""

import argparse
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle, Polygon
import numpy as np
import pandas as pd

ANALYSIS = pathlib.Path(__file__).resolve().parents[1]
WORK = ANALYSIS / "_work"
OUT = WORK / "figures" / "figure5"

# palette shared with fig1_ontology_stats.py / fig3_panels.py / notebook section 4
C_PRE, C_POST = "#8fa0a9", "#4e7e96"      # pre-workshop blue-grey, post-workshop BIH navy
C_INK, C_MUTED, C_RULE = "#003754", "#5f7078", "#c4c9ce"
C_NEWTERM, C_REPARENT = "#4e7e96", "#c8a870"   # Figure 1: new terms navy, re-parented gold
C_TINT = "#dde4e8"

UNIT_ORDER = ["clinic_letters", "icd10", "histology", "imaging", "combined"]
UNIT_LABEL = {"clinic_letters": "Clinic letters", "icd10": "ICD-10 codes",
              "histology": "Histology reports", "imaging": "Imaging reports",
              "combined": "Histology + imaging"}


def style():
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.size": 10, "font.family": "DejaVu Sans", "text.color": C_INK,
        "axes.titlesize": 11, "axes.titlepad": 8, "axes.titleweight": "normal",
        "axes.labelsize": 9.5, "axes.labelcolor": C_MUTED,
        "axes.edgecolor": C_RULE, "axes.linewidth": 1.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False,
        "xtick.color": C_MUTED, "ytick.color": C_INK, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.frameon": False, "legend.fontsize": 9,
    })


# ----------------------------------------------------------------------------- icons
def draw_icon(ax, kind, x, y, s=1.0, color=C_MUTED):
    """Minimal line icon for a data source, drawn in axes coordinates (x, y = centre)."""
    lw = 1.1
    kw = dict(transform=ax.transAxes, clip_on=False, zorder=5)
    if kind in ("clinic_letters", "icd10"):
        # a page with a folded corner; ICD-10 gets short code lines, letters get text lines
        w, h = 0.55 * s, 0.72 * s
        page = Polygon([(x - w/2, y - h/2), (x + w/2 - 0.16*s, y - h/2), (x + w/2, y - h/2 + 0.16*s),
                        (x + w/2, y + h/2), (x - w/2, y + h/2)], closed=True,
                       fc="white", ec=color, lw=lw, **kw)
        ax.add_patch(page)
        ax.plot([x + w/2 - 0.16*s, x + w/2 - 0.16*s, x + w/2], [y - h/2, y - h/2 + 0.16*s, y - h/2 + 0.16*s],
                color=color, lw=lw, **kw)
        if kind == "clinic_letters":
            for i, frac in enumerate((0.8, 0.6, 0.8, 0.45)):
                yy = y + h/2 - 0.16*s - i * 0.13*s
                ax.plot([x - w/2 + 0.1*s, x - w/2 + 0.1*s + frac * (w - 0.2*s)], [yy, yy], color=color, lw=lw, **kw)
        else:
            for i in range(3):
                yy = y + h/2 - 0.18*s - i * 0.17*s
                ax.plot([x - w/2 + 0.1*s, x - w/2 + 0.22*s], [yy, yy], color=color, lw=lw + 0.6, **kw)
                ax.plot([x - w/2 + 0.3*s, x + w/2 - 0.12*s], [yy, yy], color=color, lw=lw, **kw)
    if kind in ("histology", "combined"):
        # microscope slide with a tissue dot
        dx = -0.22 * s if kind == "combined" else 0.0
        w, h = 0.78 * s, 0.34 * s
        ax.add_patch(Rectangle((x - w/2 + dx, y - h/2), w, h, fc="white", ec=color, lw=lw, **kw))
        ax.add_patch(Rectangle((x - w/2 + dx, y - h/2), 0.18 * s, h, fc=C_TINT, ec=color, lw=lw, **kw))
        ax.add_patch(Circle((x + 0.12*s + dx, y), 0.09 * s, fc=color, ec="none", **kw))
    if kind in ("imaging", "combined"):
        # film frame with a soft blob
        dx = 0.26 * s if kind == "combined" else 0.0
        w, h = 0.6 * s, 0.6 * s
        ax.add_patch(FancyBboxPatch((x - w/2 + dx, y - h/2), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                    fc="white", ec=color, lw=lw, **kw))
        ax.add_patch(Circle((x + dx - 0.08*s, y + 0.02*s), 0.14 * s, fc=C_TINT, ec=color, lw=lw, **kw))
        ax.add_patch(Circle((x + dx + 0.1*s, y - 0.08*s), 0.1 * s, fc=C_TINT, ec=color, lw=lw, **kw))


def icon_axes(fig, bbox):
    ax = fig.add_axes(bbox)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    ax.set_aspect("equal")
    return ax


# ----------------------------------------------------------------------------- data
def load(summary_path):
    summ = pd.read_csv(summary_path).set_index("unit").loc[UNIT_ORDER]
    pat = pd.read_csv(WORK / "ehr_ic_patient_level.csv")
    terms = pd.read_csv(WORK / "ehr_ic_term_level.csv")
    return summ, pat, terms


def fmt_p(p):
    if p >= 0.001:
        return f"p = {p:.3f}"
    mant, exp = f"{p:.1e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"p = {mant} × 10{str(int(exp)).translate(sup)}"



# ----------------------------------------------------------------------------- deltas
def delta_ci(pat, unit, key, n_boot=4000, seed=20260813):
    """Mean per-patient change (post - pre) with a 95 % bootstrap percentile CI.
    mean_ic: patients scorable in both arms; sum_ic / n_terms: all patients (empty set = 0)."""
    d = pat[pat.unit == unit]
    if key == "mean_ic":
        d = d.dropna(subset=["mean_ic_pre", "mean_ic_post"])
    delta = (d[f"{key}_post"] - d[f"{key}_pre"]).values
    rng = np.random.default_rng(seed)
    boots = rng.choice(delta, (n_boot, len(delta)), replace=True).mean(axis=1)
    return delta.mean(), np.percentile(boots, 2.5), np.percentile(boots, 97.5)


def delta_panel(ax, summ, pat, key, title, fmt, ys):
    """Horizontal bars of the mean per-patient change with 95 % CI, one per data source."""
    for y, u in zip(ys, UNIT_ORDER):
        m, lo, hi = delta_ci(pat, u, key)
        ax.barh(y, m, height=.52, color=C_POST, edgecolor="white", lw=.6, zorder=2)
        ax.plot([lo, hi], [y, y], color=C_INK, lw=1.2, zorder=3)
        ax.plot([lo, lo], [y - .12, y + .12], color=C_INK, lw=1.2, zorder=3)
        ax.plot([hi, hi], [y - .12, y + .12], color=C_INK, lw=1.2, zorder=3)
        pre, post = summ.loc[u, f"{key if key != 'n_terms' else 'terms_per_patient'}_pre"], \
                    summ.loc[u, f"{key if key != 'n_terms' else 'terms_per_patient'}_post"]
        ax.text(hi, y, f"  {fmt.format(pre)} \u2192 {fmt.format(post)}", va="center", ha="left",
                fontsize=8, color=C_MUTED, zorder=4)
    ax.axvline(0, color=C_RULE, lw=1)
    ax.set_title(title, loc="left", color=C_INK)
    ax.set_yticks(ys); ax.set_ylim(-0.7, len(ys) - 0.3)
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    xmax = ax.get_xlim()[1]
    ax.set_xlim(min(0, ax.get_xlim()[0]), xmax * 1.55)
    ax.set_xlabel("\u0394 per patient (mean, 95 % CI)")

# ----------------------------------------------------------------------------- Figure 5
def figure5(summ, pat, icons=True):
    n = len(UNIT_ORDER)
    fig = plt.figure(figsize=(12, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], left=0.2 if icons else 0.14, right=0.98,
                          top=0.93, bottom=0.12, hspace=0.55, wspace=0.22)

    # a) paired boxes of per-patient mean IC per term
    ax1 = fig.add_subplot(gs[0, :])
    pos, data, cols, ticks = [], [], [], []
    for i, u in enumerate(UNIT_ORDER):
        d = pat[pat.unit == u].dropna(subset=["mean_ic_pre", "mean_ic_post"])
        data += [d.mean_ic_pre.values, d.mean_ic_post.values]
        pos += [i * 3 + 1, i * 3 + 1.85]; cols += [C_PRE, C_POST]; ticks.append(i * 3 + 1.42)
    bp = ax1.boxplot(data, positions=pos, widths=.7, patch_artist=True, showfliers=False)
    for b, c in zip(bp["boxes"], cols):
        b.set(facecolor=c, alpha=.55, edgecolor="#7c878e", lw=1.0)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color="#7c878e", lw=1.0)
    for m in bp["medians"]:
        m.set(color=C_INK, lw=1.6)
    for i, u in enumerate(UNIT_ORDER):
        ax1.text(ticks[i], 1.01, f"\u0394 = +{summ.loc[u, 'mean_ic_delta']:.3f}, N = {int(summ.loc[u, 'N']):,}\n"
                                 f"{fmt_p(summ.loc[u, 'p_mean_ic'])}",
                 transform=ax1.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.8, color=C_MUTED)
    ax1.set_xticks(ticks); ax1.set_xticklabels([UNIT_LABEL[u] for u in UNIT_ORDER], fontsize=9)
    ax1.tick_params(axis="x", length=0)
    ax1.set_ylabel("Mean IC per term, per patient")
    fig.text(0.03, ax1.get_position().y1 + 0.03, "a)", fontsize=13, fontweight="bold", color=C_INK)
    ax1.legend(handles=[mpl.patches.Patch(facecolor=C_PRE, alpha=.55, edgecolor="#7c878e", label="pre-workshop release (v2024-08-13)"),
                        mpl.patches.Patch(facecolor=C_POST, alpha=.55, edgecolor="#7c878e", label="post-workshop release (v2026-06-23)")],
               loc="upper center", bbox_to_anchor=(.5, -.2), ncol=2)
    if icons:
        y0 = ax1.get_position().y0
        for t, u in zip(ticks, UNIT_ORDER):
            fx = ax1.transData.transform((t, 0))[0] / fig.bbox.width
            draw_icon(icon_axes(fig, [fx - 0.02, y0 - 0.075, 0.04, 0.045]), u, 0.5, 0.5, s=0.9)

    # b) change per patient, three measures
    measures = [("n_terms", "Immune-branch terms per patient", "{:.2f}"),
                ("mean_ic", "Mean IC per term, per patient", "{:.3f}"),
                ("sum_ic", "Summed IC per patient", "{:.2f}")]
    ys = np.arange(n)[::-1]
    axes = [fig.add_subplot(gs[1, j]) for j in range(3)]
    for ax, (key, title, fmt) in zip(axes, measures):
        delta_panel(ax, summ, pat, key, title, fmt, ys)
        ax.set_title(title, loc="left", color=C_INK, fontsize=10)
    axes[0].set_yticklabels([UNIT_LABEL[u] for u in UNIT_ORDER], fontsize=9)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    fig.text(0.03, axes[0].get_position().y1 + 0.03, "b)", fontsize=13, fontweight="bold", color=C_INK)
    if icons:
        for y, u in zip(ys, UNIT_ORDER):
            fy = axes[0].transData.transform((0, y))[1] / fig.bbox.height
            draw_icon(icon_axes(fig, [0.03, fy - 0.03, 0.05, 0.06]), u, 0.5, 0.5, s=0.9)
    fig.text(0.5, 0.02, "b) bars: mean change per patient between the two releases (same text, same IC table); "
             "whiskers: 95 % bootstrap CI; grey: pre \u2192 post means.",
             ha="center", fontsize=8.5, color=C_MUTED)
    return fig


# ----------------------------------------------------------------------------- Supplementary Figure
def figure_s_terms(terms, icons=True, top=6):
    """Which terms were gained: post-workshop-only terms per data source."""
    units = [u for u in UNIT_ORDER if u != "combined"]   # combined would repeat histology + imaging
    n = len(units)
    fig = plt.figure(figsize=(9, 2.3 * n + 1.2))
    gs = fig.add_gridspec(n, 1, left=0.42 if icons else 0.38, right=0.97, top=0.94, bottom=0.1, hspace=0.75)
    for i, u in enumerate(units):
        ax = fig.add_subplot(gs[i, 0])
        t = terms[(terms.unit == u) & (terms.patients_pre == 0) & (terms.patients_post > 0)]
        t = t.sort_values("patients_post", ascending=False).head(top)[::-1]
        new = t.term.str.startswith("HP:521")  # ESID4HPO identifier range; others re-parented or relabelled
        ax.barh(t.label, t.patients_post, color=[C_NEWTERM if m else C_REPARENT for m in new],
                height=.62, edgecolor="white", lw=.8)
        for j, (v, ic) in enumerate(zip(t.patients_post, t.ic)):
            ax.text(v + t.patients_post.max() * 0.02, j, f"{int(v)}  (IC {ic:.2f})", va="center", fontsize=8, color=C_MUTED)
        ax.set_xlim(0, t.patients_post.max() * 1.45)
        ax.set_title(f"{UNIT_LABEL[u]}: terms recoverable only with the post-workshop release",
                     loc="left", color=C_INK, fontsize=10)
        ax.tick_params(axis="y", length=0, labelsize=8.5); ax.spines["left"].set_visible(False)
        ax.set_xlabel("patients" if i == n - 1 else "")
        if icons:
            pos = ax.get_position()
            draw_icon(icon_axes(fig, [0.02, pos.y0 + pos.height / 2 - 0.035, 0.06, 0.07]), u, 0.5, 0.5, s=0.9)
    handles = [mpl.patches.Patch(facecolor=C_NEWTERM, label="term created by ESID4HPO"),
               mpl.patches.Patch(facecolor=C_REPARENT, label="pre-existing term, re-parented into the immune branch or relabelled")]
    fig.legend(handles=handles, loc="lower center", ncol=1, bbox_to_anchor=(0.6, 0.0))
    return fig


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=pathlib.Path, default=WORK / "ehr_ic_summary_sa.csv",
                    help="summary table to take N, means and p from (default: Alex's numbers)")
    ap.add_argument("--only", choices=["main", "supp"], help="draw only Figure 5 or only the supplementary figure")
    ap.add_argument("--no-icons", action="store_true")
    args = ap.parse_args()
    style()
    summ, pat, terms = load(args.summary)
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = {"main": ("fig5_ehr", lambda: figure5(summ, pat, not args.no_icons)),
            "supp": ("figS_ehr_terms", lambda: figure_s_terms(terms, not args.no_icons))}
    for key, (name, build) in jobs.items():
        if args.only and key != args.only:
            continue
        fig = build()
        for ext in ("png", "pdf", "svg"):
            fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUT / name}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
