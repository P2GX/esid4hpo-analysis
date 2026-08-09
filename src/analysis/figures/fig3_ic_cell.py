"""
Drop-in replacement for the intrinsic-IC cell in esid4hpo_analysis.ipynb.

Method (revised per Peter's comment "HPOAs durch Kohorte ersetzen, und dann in
Kohorte den IC der Terms bestimmen"): Resnik IC, IC(t) = -ln p(t), with term
frequencies taken from THE COHORT ITSELF rather than from phenotype.hpoa.

    freq(t) = number of cohort patients annotated to t OR to any descendant of t
    p(t)    = freq(t) / n_patients

Output schema is unchanged, so fig3_panels.py needs no edit:

    semantic_term_level.csv     cohort, phenopacket, version, curie, ic
    semantic_patient_level.csv  cohort, phenopacket, n_new, n_aged, ic_new,
                                ic_aged, mean_ic_new, mean_ic_aged,
                                n_newly_enabled, frac_newly_enabled

Expected inputs already in the notebook
---------------------------------------
    hpo_new     MinimalOntology, v2026-06-23      (fixed scoring graph)
    terms_new   {cohort: {phenopacket: [curie, ...]}}   curated arm
    terms_aged  {cohort: {phenopacket: [curie, ...]}}   aged arm
Rename on the lines marked ADJUST if your variables differ.
"""

# ---------------------------------------------------------------- cell start
import pandas as pd
from ic_annotation import AncestorIndex, build_ic_model, score_annotation_set

# --- configuration ---------------------------------------------------------
# "cohort" : frequencies from that cohort's patients only.   <-- PRIMARY
#            This is Peter's specification. n is small (21-38), so IC is
#            coarse and MUST be normalised for cohorts to be comparable.
# "pooled" : frequencies from all three cohorts' patients pooled (n = 87).
#            Still only "what we actually annotated" - no HPOA involved.
#            Report as a sensitivity analysis; see the note at the bottom.
IC_CORPUS    = "cohort"
IC_NORMALIZE = True    # divide by ln(n) -> IC on [0,1]; required for "cohort"
IC_SMOOTHING = 0.0     # not needed: corpus covers both arms by construction

ONTO  = hpo_new                                   # ADJUST: fixed graph v2026-06-23
TERMS = {"new": terms_new, "aged": terms_aged}    # ADJUST

index = AncestorIndex(ONTO)
COHORTS_IN_SCOPE = list(TERMS["new"].keys())


def _entities(cohorts):
    """One entity per patient; union of both arms.

    Using the union guarantees every scored term has frequency >= 1 so no
    smoothing is needed. Under pure ancestor-ageing the union is identical to
    the curated arm alone (aged terms are ancestors, already in the closure),
    so this costs nothing and is robust to dropped excluded features.
    """
    ents = []
    for coh in cohorts:
        for pp in TERMS["new"][coh]:
            ents.append(set(TERMS["new"][coh][pp]) | set(TERMS["aged"][coh].get(pp, [])))
    return ents


# --- build ONE frozen frequency table per scoring unit ---------------------
# Both arms are scored with the same table. Building a table per arm would put
# the arms in different probability spaces and the IC difference would no
# longer be interpretable. This is the direct analogue of the fixed-graph
# argument and is the single most important thing to keep right.
if IC_CORPUS == "cohort":
    models = {
        coh: build_ic_model(_entities([coh]), index, corpus="cohort",
                            normalize=IC_NORMALIZE, smoothing=IC_SMOOTHING)
        for coh in COHORTS_IN_SCOPE
    }
elif IC_CORPUS == "pooled":
    m = build_ic_model(_entities(COHORTS_IN_SCOPE), index, corpus="pooled",
                       normalize=IC_NORMALIZE, smoothing=IC_SMOOTHING)
    models = {coh: m for coh in COHORTS_IN_SCOPE}
else:
    raise ValueError(IC_CORPUS)

for coh, m in models.items():
    print(coh, m.describe())
    if m.unscoreable:
        print("   not under HP:0000118 / unknown to graph:", sorted(m.unscoreable)[:8])

# --- score -----------------------------------------------------------------
term_rows, pat_rows = [], []
for coh in COHORTS_IN_SCOPE:
    m = models[coh]
    for pp, t_new in TERMS["new"][coh].items():
        t_aged = TERMS["aged"][coh].get(pp, [])

        s_new,  mu_new,  n_new,  _ = score_annotation_set(t_new,  m, index)
        s_aged, mu_aged, n_aged, _ = score_annotation_set(t_aged, m, index)

        for version, tt in (("new", t_new), ("aged", t_aged)):
            for c in tt:
                pid = index.primary_id(c)
                if pid in m.ic:
                    term_rows.append({"cohort": coh, "phenopacket": pp,
                                      "version": version, "curie": pid,
                                      "ic": m.ic[pid]})

        newly = len(set(t_new) - set(t_aged))
        pat_rows.append({
            "cohort": coh, "phenopacket": pp,
            "n_new": n_new, "n_aged": n_aged,
            "ic_new": s_new, "ic_aged": s_aged,
            "mean_ic_new": mu_new, "mean_ic_aged": mu_aged,
            "n_newly_enabled": newly,
            "frac_newly_enabled": newly / len(t_new) if t_new else float("nan"),
        })

semantic_term_level    = pd.DataFrame(term_rows)
semantic_patient_level = pd.DataFrame(pat_rows)
semantic_term_level.to_csv(FPATH_WORK / "semantic_term_level.csv", index=False)
semantic_patient_level.to_csv(FPATH_WORK / "semantic_patient_level.csv", index=False)

# --- diagnostic: which new terms are invisible to a cohort corpus? ---------
# A term carried by EVERY patient has p = 1 and therefore IC = 0, in both arms.
# So a workshop term that is pathognomonic for the cohort contributes ZERO
# measured gain under corpus="cohort". Print them so the effect is explicit
# rather than silent - this is the main caveat of Peter's specification.
for coh in COHORTS_IN_SCOPE:
    m = models[coh]
    saturated = sorted(t for t, f in m.freq.items() if f == m.n_entities)
    direct = {index.primary_id(c) for pp in TERMS["new"][coh]
              for c in TERMS["new"][coh][pp]}
    hit = [t for t in saturated if t in direct]
    if hit:
        print(f"{coh}: {len(hit)} directly-annotated terms carried by all "
              f"{m.n_entities} patients -> IC = 0, no measurable gain")
        for t in hit[:10]:
            lbl = ONTO.get_term(t)
            print("   ", t, lbl.name if lbl else "")

semantic_patient_level.head()
# ------------------------------------------------------------------ cell end
