# =====================================================================
# FIGURE 3 - information content (panels A-D)
#   A  paired per-patient IC per cohort (+ Holm-corrected p-values)
#   B  cross-cohort summary bars with % gain
#   C  fraction of annotations using post-workshop terms
#   D  worked example: annotation collapse in a single patient
#
# Panels are saved separately (png/pdf/svg) in FPATH_FIGS for assembly
# in PowerPoint. Cohorts listed in FIG3_COHORTS but not yet analysed
# (e.g. NFKB1) still get a reserved, labelled slot so the figure layout
# does not change when their data arrive.
#
# Run from the notebook with:
#     exec(open(FIG3_PATH).read())
# (the loader cell resolves FIG3_PATH robustly - see the notebook)
# =====================================================================
try:
    from scipy import stats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False
    print("scipy unavailable - p-values will be skipped")

# cohorts to show, in order; those without data get a reserved empty slot
FIG3_COHORTS = ["SOCS1", "APDS", "NFKB1"]

# which panels to produce: "A" paired IC, "B" summary bars, "C" post-workshop
# fraction, "D" annotation collapse. Panel B is redundant with A (same data,
# and the % gain is now printed on A) so it is off by default - add "B" back
# if you want it as a separate panel.
FIG3_PANELS = ["A", "C", "D"]

# Panels A and C are stacked in the final figure, so they are rendered on a
# fixed canvas of identical width. NOTE: they are saved WITHOUT
# bbox_inches="tight" - tight cropping trims each figure to its own content and
# would make the widths differ again. Layout is handled by tight_layout()
# inside the canvas instead.
FIG3_WIDTH = 9.6   # inches, shared by panels A and C
FIG3_H_A   = 4.6   # panel A height
FIG3_H_C   = 2.9   # panel C height (flatter)

C_OLDBAR = "#dde4e8"   # pre-workshop  (BIH navy tint)
C_NEWBAR = "#4e7e96"   # post-workshop (BIH navy)
C_ACC    = "#c8a870"   # accent: merged / collapsed
C_INK    = "#003754"
C_MUTED  = "#5f7078"
C_PEND   = "#aab4ba"   # pending / reserved slot


def _f3style():
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": None,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.size": 10.5, "font.family": "DejaVu Sans", "text.color": C_INK,
        "axes.titlesize": 11.5, "axes.titlepad": 8,
        "axes.labelsize": 10, "axes.labelcolor": C_MUTED,
        "axes.edgecolor": "#c4c9ce", "axes.linewidth": 1.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": C_MUTED, "ytick.color": C_INK,
        "legend.frameon": False, "legend.fontsize": 9,
    })


def save_panel(fig, name, tight=False):
    """tight=False keeps the exact figsize (so stacked panels share a width);
    tight=True crops to content (used for the wide panel D)."""
    FPATH_FIGS.mkdir(parents=True, exist_ok=True)
    kw = dict(bbox_inches="tight") if tight else {}
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FPATH_FIGS / f"{name}.{ext}", **kw)
    w, h = fig.get_size_inches()
    print(f"  saved {name}.png/.pdf/.svg  ({w:.2f} x {h:.2f} in)")


def _cohort_data(df_patient, ch):
    """Rows for a cohort; empty DataFrame if it has no data yet."""
    return df_patient[df_patient.cohort == ch]


def _pending(ax, ch):
    """Draw a reserved, clearly-labelled empty slot for a cohort with no data."""
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linestyle((0, (4, 4)))
        sp.set_color("#d3d9dd")
    ax.text(0.5, 0.5, f"{ch}\n(pending)", ha="center", va="center",
            fontsize=11, color=C_PEND, transform=ax.transAxes)


# ---------------------------------------------------------------- stats
def _rank_biserial(new_v, old_v):
    d = np.asarray(new_v, float) - np.asarray(old_v, float)
    d = d[d != 0]
    if d.size == 0:
        return np.nan
    r = stats.rankdata(np.abs(d))
    return (r[d > 0].sum() - r[d < 0].sum()) / r.sum()


def _holm(p):
    p = np.asarray(p, float)
    m = p.size
    adj = np.empty(m)
    run = 0.0
    for k, i in enumerate(np.argsort(p)):
        run = max(run, (m - k) * p[i])
        adj[i] = min(run, 1.0)
    return adj


def ic_stats(df_patient, alternative="two-sided"):
    """Paired Wilcoxon per cohort + effect size + Holm correction across cohorts."""
    rows = []
    for ch in FIG3_COHORTS:
        s = _cohort_data(df_patient, ch).dropna(subset=["mean_ic_aged", "mean_ic_new"])
        if not len(s):
            rows.append(dict(cohort=ch, n=0, n_informative=0, W=np.nan,
                             p_raw=np.nan, rank_biserial=np.nan, gain_pct=np.nan))
            continue
        aged, newv = s.mean_ic_aged.values, s.mean_ic_new.values
        n_inf = int(np.sum(newv - aged != 0))
        if _HAVE_SCIPY and n_inf:
            W, p = stats.wilcoxon(newv, aged, alternative=alternative)
        else:
            W, p = np.nan, np.nan
        rows.append(dict(cohort=ch, n=len(s), n_informative=n_inf, W=W, p_raw=p,
                         rank_biserial=_rank_biserial(newv, aged) if _HAVE_SCIPY else np.nan,
                         gain_pct=((newv.mean() - aged.mean()) / aged.mean() * 100
                                   if aged.mean() else np.nan)))
    res = pd.DataFrame(rows)
    m = res.p_raw.notna()
    res["p_holm"] = np.nan
    if m.any():
        res.loc[m, "p_holm"] = _holm(res.loc[m, "p_raw"].values)
    return res


def _stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n.s."
    return ("****" if p < 1e-4 else "***" if p < 1e-3
            else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s.")


# ------------------------------------------------------------- panel A
def fig3_panelA(df_patient, stats_tbl):
    """Paired per-patient IC. Clean strip + box: no connecting lines, points sit
    exactly on the category centre (overlap shown by transparency)."""
    _f3style()
    n = len(FIG3_COHORTS)
    fig, axes = plt.subplots(1, n, figsize=(FIG3_WIDTH, FIG3_H_A),
                             squeeze=False, sharey=True)
    for ax, ch in zip(axes[0], FIG3_COHORTS):
        s = _cohort_data(df_patient, ch).dropna(subset=["mean_ic_aged", "mean_ic_new"])
        if not len(s):
            _pending(ax, ch)
            ax.set_title(ch, fontweight="bold", color=C_PEND)
            continue

        bp = ax.boxplot([s.mean_ic_aged, s.mean_ic_new], positions=[1, 2],
                        widths=.46, showfliers=False, patch_artist=True, zorder=2)
        for b_, c_ in zip(bp["boxes"], [C_OLDBAR, C_NEWBAR]):
            b_.set(facecolor=c_, alpha=.55, edgecolor="#7c878e", lw=1.0)
        for el in bp["whiskers"] + bp["caps"]:
            el.set(color="#7c878e", lw=1.0)
        for el in bp["medians"]:
            el.set(color=C_INK, lw=1.8)

        # points exactly on the category centre (no jitter, no linking lines)
        ax.scatter(np.ones(len(s)), s.mean_ic_aged, s=22, color="#8fa0a9",
                   alpha=.55, linewidths=0, zorder=3)
        ax.scatter(np.full(len(s), 2.0), s.mean_ic_new, s=22, color=C_NEWBAR,
                   alpha=.55, linewidths=0, zorder=3)

        r = stats_tbl.loc[stats_tbl.cohort == ch]
        if len(r) and not np.isnan(r.iloc[0].p_holm):
            r = r.iloc[0]
            lo = min(s.mean_ic_aged.min(), s.mean_ic_new.min())
            hi = max(s.mean_ic_aged.max(), s.mean_ic_new.max())
            pad = (hi - lo) * 0.12 if hi > lo else 0.05
            yb = hi + pad * 0.55
            ax.plot([1, 1, 2, 2], [yb, yb + pad * .28, yb + pad * .28, yb],
                    color=C_MUTED, lw=1.0, zorder=4)
            gain = "" if np.isnan(r.gain_pct) else f"   (+{r.gain_pct:.1f}%)"
            ax.text(1.5, yb + pad * .34,
                    f"{_stars(r.p_holm)}  p = {r.p_holm:.1e}{gain}",
                    ha="center", va="bottom", fontsize=8.5, color=C_INK)
            ax.set_ylim(lo - pad * .6, hi + pad * 1.5)

        ax.set_xlim(0.45, 2.55)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"pre\n{HPO_OLD_TAG}", f"post\n{HPO_NEW_TAG}"])
        ax.set_title(f"{ch} (n = {len(s)})", fontweight="bold", color=C_INK)
        ax.grid(axis="y", color="#000", alpha=.06, lw=.8)
        ax.set_axisbelow(True)
        ax.tick_params(length=0)
    axes[0][0].set_ylabel("Mean information content per term")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------- panel B
def fig3_panelB(df_patient, stats_tbl):
    _f3style()
    xs = np.arange(len(FIG3_COHORTS))
    w = .36
    fig, ax = plt.subplots(figsize=(1.9 * len(FIG3_COHORTS) + 2.4, 4.5))
    for off, col, lab, key in [
            (-w / 2, C_OLDBAR, f"pre-workshop ({HPO_OLD_TAG})", "mean_ic_aged"),
            (w / 2, C_NEWBAR, f"post-workshop ({HPO_NEW_TAG})", "mean_ic_new")]:
        means, errs = [], []
        for ch in FIG3_COHORTS:
            sub = _cohort_data(df_patient, ch)
            v = sub[key].dropna() if len(sub) else pd.Series(dtype=float)
            means.append(v.mean() if len(v) else 0.0)
            errs.append(1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)
        ax.bar(xs + off, means, w, yerr=errs, capsize=3, color=col,
               edgecolor="white", lw=.8, label=lab,
               error_kw=dict(lw=1, ecolor="#6b7279"))
    ymax = max(list(ax.get_ylim()) + [0.01])
    for i, ch in enumerate(FIG3_COHORTS):
        s = _cohort_data(df_patient, ch)
        r = stats_tbl.loc[stats_tbl.cohort == ch]
        if not len(s):
            ax.text(i, ymax * 0.06, f"{ch}\n(pending)", ha="center", va="bottom",
                    fontsize=9.5, color=C_PEND)
            continue
        if len(r) and not np.isnan(r.iloc[0].gain_pct):
            top = max(s.mean_ic_new.mean(), s.mean_ic_aged.mean())
            ax.text(i, top * 1.12, f"+{r.iloc[0].gain_pct:.0f}%", ha="center",
                    fontsize=10.5, fontweight="bold", color=C_NEWBAR)
    ax.set_xticks(xs)
    ax.set_xticklabels(FIG3_COHORTS)
    ax.set_ylabel("Mean information content per term")
    ax.grid(axis="y", color="#000", alpha=.06, lw=.8)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.09), ncol=2)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------- panel C
def fig3_panelC(df_patient):
    _f3style()
    fig, ax = plt.subplots(figsize=(FIG3_WIDTH, FIG3_H_C))
    vals, errs, have = [], [], []
    for ch in FIG3_COHORTS:
        s = _cohort_data(df_patient, ch)
        v = s.frac_newly_enabled.dropna() * 100 if len(s) else pd.Series(dtype=float)
        have.append(bool(len(v)))
        vals.append(v.mean() if len(v) else 0.0)
        errs.append(1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)
    bars = ax.bar(FIG3_COHORTS, vals, .42, yerr=errs, capsize=3, color=C_NEWBAR,
                  edgecolor="white", lw=.8, error_kw=dict(lw=1, ecolor="#6b7279"))
    top = max(vals) if any(have) else 1.0
    for b, v, e, h, ch in zip(bars, vals, errs, have, FIG3_COHORTS):
        if h:
            # value label INSIDE the bar, near its base. It cannot sit just under
            # the bar top because the error bar reaches down into the bar there
            # and would cross the text.
            if v > top * 0.20:
                ax.text(b.get_x() + b.get_width() / 2, top * .06, f"{v:.1f}%",
                        ha="center", va="bottom", fontsize=10.5, fontweight="bold",
                        color="white")
            else:
                # bar too short for inside text -> above the upper whisker cap
                ax.text(b.get_x() + b.get_width() / 2, v + e + top * .05,
                        f"{v:.1f}%", ha="center", va="bottom", fontsize=10.5,
                        fontweight="bold", color=C_INK)
        else:
            ax.text(b.get_x() + b.get_width() / 2, top * .06, f"{ch}\n(pending)",
                    ha="center", va="bottom", fontsize=9.5, color=C_PEND)
    # short label: a long rotated ylabel is taller than the flat axes and clips
    ax.set_ylabel("Post-workshop terms (%)")
    ax.set_ylim(0, top * 1.25)
    ax.margins(x=0.06)
    ax.grid(axis="y", color="#000", alpha=.06, lw=.8)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.subplots_adjust(top=0.94, bottom=0.16)
    return fig


# ------------------------------------------------------------- panel D
def collapse_map(pp, hpo_old, hpo_new, include_excluded=True):
    return [(c, e, age_feature(c, e, hpo_old, hpo_new))
            for c, e in feature_pairs(pp, include_excluded)]


def pick_collapse_example(cohort, hpo_old, hpo_new, include_excluded=True):
    best = None
    for fname, pp in load_phenopackets(cohort.pp_dir()):
        m = collapse_map(pp, hpo_old, hpo_new, include_excluded)
        groups = {}
        for c, e, a in m:
            if a is not None:
                groups.setdefault(a, []).append(c)
        merged = sum(len(v) - 1 for v in groups.values() if len(v) > 1)
        dropped = sum(1 for _, _, a in m if a is None)
        score = merged * 2 + dropped
        if best is None or score > best[0]:
            best = (score, fname, m, groups)
    return best


def _short(txt, n=46):
    return txt if len(txt) <= n else txt[:n - 1] + "\u2026"


def fig3_panelD(cohort, hpo_old, hpo_new, max_rows=14):
    _f3style()
    picked = pick_collapse_example(cohort, hpo_old, hpo_new)
    if picked is None:
        print("  no phenopackets for", cohort.name)
        return None, None
    score, fname, m, groups = picked
    rows = [(c, e, a) for c, e, a in m if (a is None) or (a != c)]
    rows.sort(key=lambda t: (t[2] is None, t[2] or ""))
    rows = rows[:max_rows]
    if not rows:
        print("  no collapse to illustrate in", cohort.name)
        return None, fname

    src = {}
    for i, (_, _, a) in enumerate(rows):
        src.setdefault(a if a is not None else "__dropped__", []).append(i)

    fig, ax = plt.subplots(figsize=(11.6, 0.52 * len(rows) + 1.7))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    yL = {i: 1 - (i + 0.5) / len(rows) for i in range(len(rows))}
    yR = {k: float(np.mean([yL[i] for i in v])) for k, v in src.items()}
    XL, WL, XR, WR = 0.005, 0.435, 0.560, 0.435

    def box(x, y, w, text, fc, ec, tc, bold=False):
        ax.add_patch(mpl.patches.FancyBboxPatch(
            (x, y - 0.026), w, 0.052,
            boxstyle="round,pad=0.004,rounding_size=0.010",
            facecolor=fc, edgecolor=ec, lw=1.2))
        ax.text(x + 0.012, y, text, va="center", ha="left", fontsize=8.5,
                color=tc, fontweight="bold" if bold else "normal")

    for i, (c, e, a) in enumerate(rows):
        key = a if a is not None else "__dropped__"
        merged = len(src[key]) > 1
        lab = _short(hpo_new.get_term_name(c))
        if e:
            lab = "excluded: " + lab
        box(XL, yL[i], WL, lab, C_NEWBAR, C_NEWBAR, "white")
        ax.annotate("", xy=(XR - 0.004, yR[key]), xytext=(XL + WL + 0.004, yL[i]),
                    arrowprops=dict(arrowstyle="-|>",
                                    color=C_ACC if merged else "#b9c3c9",
                                    lw=1.7 if merged else 1.0,
                                    shrinkA=0, shrinkB=0))
    for key, ys in yR.items():
        n = len(src[key])
        if key == "__dropped__":
            box(XR, ys, WR, "not representable \u2014 annotation lost",
                "white", "#c4c9ce", C_MUTED)
        else:
            lab = _short(hpo_old.get_term_name(key))
            suffix = f"    [{n} terms merged]" if n > 1 else ""
            box(XR, ys, WR, lab + suffix, C_OLDBAR,
                C_ACC if n > 1 else "#b9c3c9", C_INK, bold=n > 1)

    ax.text(XL, 1.015, f"Curated with {HPO_NEW_TAG}", fontsize=10.5,
            fontweight="bold", color=C_INK, va="bottom")
    ax.text(XR, 1.015, f"Same patient expressed in {HPO_OLD_TAG}", fontsize=10.5,
            fontweight="bold", color=C_INK, va="bottom")
    ax.text(0, -0.035, f"{cohort.name}  \u00b7  {fname}", fontsize=9, color=C_MUTED)
    fig.tight_layout()
    return fig, fname


# ------------------------------------------------------------- driver
print(f"Figure 3 -> {FPATH_FIGS}")
stats_tbl = ic_stats(df_patient, alternative="two-sided")
try:
    display(stats_tbl.round(4))
except NameError:
    print(stats_tbl.round(4))
FPATH_WORK.mkdir(parents=True, exist_ok=True)
stats_tbl.to_csv(FPATH_WORK / "fig3_ic_stats.csv", index=False)

if "A" in FIG3_PANELS:
    save_panel(fig3_panelA(df_patient, stats_tbl), "fig3A_paired_ic")
if "B" in FIG3_PANELS:
    save_panel(fig3_panelB(df_patient, stats_tbl), "fig3B_summary_gain")
if "C" in FIG3_PANELS:
    save_panel(fig3_panelC(df_patient), "fig3C_post_workshop_fraction")

if "D" in FIG3_PANELS:
    _with_data = [c for c in ACTIVE if len(_cohort_data(df_patient, c.name))]
    if _with_data:
        figD, ex_name = fig3_panelD(_with_data[0], hpo_old, hpo_new)
        if figD is not None:
            save_panel(figD, "fig3D_annotation_collapse", tight=True)
            print("  collapse example:", _with_data[0].name, ex_name)
print("Figure 3 done.")
