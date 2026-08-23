# =====================================================================
# FIGURE 3 - information content (panels A and C)
#   A  paired per-patient IC per cohort (mean paired delta-IC with 95% CI;
#      Wilcoxon/Holm stays in the stats table, not on the panel)
#   C  fraction of annotations using post-workshop terms
# (panels B and D removed 2026-08-23; recover from git history if needed)
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
FIG3_COHORTS = ["SOCS1", "APDS1", "NFKB1"]

# which panels to produce: "A" paired IC, "C" post-workshop fraction
FIG3_PANELS = ["A", "C"]

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
    """tight=False keeps the exact figsize (so stacked panels A and C share a width)."""
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
                             p_raw=np.nan, rank_biserial=np.nan, gain_pct=np.nan,
                             delta_mean=np.nan, delta_lo=np.nan, delta_hi=np.nan,
                             hl_delta=np.nan, hl_lo=np.nan, hl_hi=np.nan,
                             n_improved=0, n_worse=0))
            continue
        aged, newv = s.mean_ic_aged.values, s.mean_ic_new.values
        n_inf = int(np.sum(newv - aged != 0))
        if _HAVE_SCIPY and n_inf:
            W, p = stats.wilcoxon(newv, aged, alternative=alternative)
        else:
            W, p = np.nan, np.nan
        d = newv - aged
        boot_lo, boot_hi = _boot_ci_mean(d)
        hl, hl_lo, hl_hi = _hodges_lehmann(d) if _HAVE_SCIPY else (np.nan,) * 3
        rows.append(dict(cohort=ch, n=len(s), n_informative=n_inf, W=W, p_raw=p,
                         rank_biserial=_rank_biserial(newv, aged) if _HAVE_SCIPY else np.nan,
                         gain_pct=((newv.mean() - aged.mean()) / aged.mean() * 100
                                   if aged.mean() else np.nan),
                         delta_mean=float(d.mean()), delta_lo=boot_lo, delta_hi=boot_hi,
                         hl_delta=hl, hl_lo=hl_lo, hl_hi=hl_hi,
                         n_improved=int((d > 0).sum()), n_worse=int((d < 0).sum())))
    res = pd.DataFrame(rows)
    m = res.p_raw.notna()
    res["p_holm"] = np.nan
    if m.any():
        res.loc[m, "p_holm"] = _holm(res.loc[m, "p_raw"].values)
    return res


def _boot_ci_mean(d, n_boot=10000, seed=20260813, level=0.95):
    """Bootstrap percentile CI for the mean paired difference. Fixed seed so
    the figure is reproducible run-to-run."""
    d = np.asarray(d, float)
    if d.size == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    means = d[idx].mean(axis=1)
    a = (1.0 - level) / 2.0
    return float(np.quantile(means, a)), float(np.quantile(means, 1.0 - a))


def _hodges_lehmann(d, level=0.95):
    """Pseudomedian (Hodges-Lehmann) of the paired differences with a CI from
    the Walsh averages, normal approximation to the signed-rank distribution
    (what R's wilcox.test(conf.int=TRUE) computes). Zero differences are
    retained, so with many unchanged patients the estimate can be 0 even when
    the Wilcoxon on the informative pairs is significant (NFKB1!) - that is
    honest, not a bug; report n_improved alongside. Written to the stats CSV
    as a companion to delta_mean."""
    d = np.asarray(d, float)
    n = d.size
    if n == 0:
        return np.nan, np.nan, np.nan
    walsh = np.sort(np.array([(d[i] + d[j]) / 2.0
                              for i in range(n) for j in range(i, n)]))
    est = float(np.median(walsh))
    m = walsh.size
    z = stats.norm.ppf(1.0 - (1.0 - level) / 2.0)
    k = int(np.floor(n * (n + 1) / 4.0
                     - z * np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)))
    k = max(k, 0)
    if k >= m - k - 1:
        return est, float(walsh[0]), float(walsh[-1])
    return est, float(walsh[k]), float(walsh[m - k - 1])


# ------------------------------------------------------------- panel A
def fig3_panelA(df_patient, stats_tbl):
    """Paired per-patient IC. Clean strip + box: no connecting lines, points sit
    exactly on the category centre (overlap shown by transparency).

    The annotation reports the mean paired difference (aged -> curated) with a
    95% bootstrap CI plus how many patients improved. The direction of the
    difference is fixed by construction (ageing can only lower IC), so no
    p-value is printed on the panel; the Wilcoxon tests remain in the stats
    table (fig3_ic_stats.csv / Table 1)."""
    _f3style()
    n = len(FIG3_COHORTS)
    fig, axes = plt.subplots(1, n, figsize=(FIG3_WIDTH, FIG3_H_A),
                             squeeze=False, sharey=True)

    # Global y-range, computed once: with sharey=True a per-axes set_ylim is
    # overridden by whichever cohort is drawn last (the old behaviour - the
    # visible limits were silently those of the final cohort).
    have = df_patient.dropna(subset=["mean_ic_aged", "mean_ic_new"])
    have = have[have.cohort.isin(FIG3_COHORTS)]
    if len(have):
        glo = float(min(have.mean_ic_aged.min(), have.mean_ic_new.min()))
        ghi = float(max(have.mean_ic_aged.max(), have.mean_ic_new.max()))
    else:
        glo, ghi = 0.0, 1.0
    gpad = (ghi - glo) * 0.12 if ghi > glo else 0.05

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
        if len(r) and not np.isnan(r.iloc[0].delta_mean):
            r = r.iloc[0]
            hi = max(s.mean_ic_aged.max(), s.mean_ic_new.max())
            yb = hi + gpad * 0.55
            ax.plot([1, 1, 2, 2], [yb, yb + gpad * .28, yb + gpad * .28, yb],
                    color=C_MUTED, lw=1.0, zorder=4)
            gain = "" if np.isnan(r.gain_pct) else f"  (+{r.gain_pct:.1f}%)"
            line1 = (f"\u0394IC = +{r.delta_mean:.3f} "
                     f"[{r.delta_lo:.3f}, {r.delta_hi:.3f}]{gain}")
            line2 = f"{int(r.n_improved)}/{int(r.n)} patients improved"
            ax.text(1.5, yb + gpad * .34, line1 + "\n" + line2,
                    ha="center", va="bottom", fontsize=8.0, color=C_INK)

        ax.set_xlim(0.45, 2.55)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"pre\n{HPO_OLD_TAG}", f"post\n{HPO_NEW_TAG}"])
        ax.set_title(f"{ch} (n = {len(s)})", fontweight="bold", color=C_INK)
        ax.grid(axis="y", color="#000", alpha=.06, lw=.8)
        ax.set_axisbelow(True)
        ax.tick_params(length=0)
    axes[0][0].set_ylim(glo - gpad * .6, ghi + gpad * 2.6)
    axes[0][0].set_ylabel("Mean information content per term")
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
if "C" in FIG3_PANELS:
    save_panel(fig3_panelC(df_patient), "fig3C_post_workshop_fraction")
print("Figure 3 done.")
