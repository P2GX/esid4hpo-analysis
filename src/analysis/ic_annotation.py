"""
ESID4HPO -- annotation-based information content (revised Figure 3 method).

Background
----------
The first draft scored terms with *intrinsic* IC (Seco et al. 2004),
IC(t) = 1 - ln(|desc(t)|)/ln(N). That measure depends only on where a term sits
in the graph, so it really reports how finely the ontology happens to be
subdivided in a given region -- not how informative the term is about a patient.
It was also insensitive by construction to re-parenting and redefinition.

This module implements *Resnik* IC instead,

    IC(t) = -ln p(t),        p(t) = freq(t) / |C|

but takes the term frequencies from the **cohort itself** rather than from
`phenotype.hpoa`. That was the point of the objection: the stock HPOA does not
yet carry the post-workshop immunophenotyping terms, so frequencies drawn from
it are undefined exactly where the analysis needs them.

`freq(t)` counts corpus entities annotated to `t` **or to any descendant of
`t`** -- i.e. each entity's direct annotations are expanded to their ancestor
closure before counting. This is the standard true-path-rule propagation.

Three design knobs matter scientifically and are therefore explicit
parameters rather than buried defaults. See `CORPUS_CHOICES` and the notes on
`build_ic_model`.

Author: ESID4HPO analysis
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

__all__ = [
    "PHENOTYPIC_ABNORMALITY",
    "CORPUS_CHOICES",
    "AncestorIndex",
    "ICModel",
    "build_ic_model",
    "score_annotation_set",
    "load_hpoa_disease_annotations",
    "build_hpoa_replaced_corpus",
]

PHENOTYPIC_ABNORMALITY = "HP:0000118"

CORPUS_CHOICES = (
    # entities = the patients of the cohort under analysis.
    # Literal reading of "replace the HPOA with the cohort". IC is bounded by
    # ln(n), so with n=28 the ceiling is 3.33 nats and resolution is coarse.
    # Normalisation (below) is essentially mandatory here, otherwise a larger
    # cohort mechanically yields larger IC and cohorts are not comparable.
    "cohort",
    # entities = patients pooled across all cohorts.
    # One common probability space for every cohort, so cross-cohort comparison
    # is meaningful and the ceiling rises to ln(n_total). Recommended default.
    "pooled",
    # entities = every disease in phenotype.hpoa, with the target disease's own
    # annotation block deleted and the cohort's patients inserted as entities.
    # Keeps the ~8k-disease background that makes Resnik IC stable, while
    # letting the new terms acquire non-zero frequency.
    "hpoa_replaced",
)


# --------------------------------------------------------------------------
# ancestor closure
# --------------------------------------------------------------------------
class AncestorIndex:
    """Memoised ancestor lookup over a fixed hpotk ontology.

    All arms of the analysis must be scored against a *single* ontology graph
    (v2026-06-23), so that the only thing varying between arms is which term the
    patient carries. This class is the one place the graph is touched.
    """

    def __init__(self, ontology, root: str = PHENOTYPIC_ABNORMALITY):
        self._onto = ontology
        self._root = root
        self._cache: Dict[str, Optional[Set[str]]] = {}
        self._primary: Dict[str, Optional[str]] = {}

    # -- id handling -------------------------------------------------------
    def primary_id(self, curie: str) -> Optional[str]:
        """Resolve a CURIE to its primary (non-obsolete, non-alt) form.

        Returns None if the term is unknown to this release -- callers decide
        whether that is an error or an expected drop.
        """
        if curie in self._primary:
            return self._primary[curie]
        term = self._onto.get_term(curie)
        pid = None if term is None else term.identifier.value
        self._primary[curie] = pid
        return pid

    # -- closure -----------------------------------------------------------
    def ancestors(self, curie: str, include_self: bool = True) -> Set[str]:
        """Ancestors of `curie` restricted to the `root` subtree.

        The root itself is included when it is an ancestor; `include_self` adds
        the term. Unknown terms yield an empty set.
        """
        key = f"{curie}|{int(include_self)}"
        cached = self._cache.get(key)
        if cached is not None:
            return cached

        pid = self.primary_id(curie)
        if pid is None:
            self._cache[key] = set()
            return set()

        tid = self._onto.get_term(pid).identifier
        anc = {a.value for a in self._onto.graph.get_ancestors(tid)}
        if include_self:
            anc.add(pid)

        # keep only the phenotypic-abnormality subtree; terms outside it
        # (clinical modifiers, frequency, inheritance) are not scored
        if self._root is not None:
            if self._root not in anc and pid != self._root:
                # term is not under the root at all -> not scoreable
                anc = set()
            else:
                anc = {a for a in anc if a == self._root or self._is_under_root(a)}

        self._cache[key] = anc
        return anc

    def _is_under_root(self, curie: str) -> bool:
        if curie == self._root:
            return True
        tid = self._onto.get_term(curie)
        if tid is None:
            return False
        return self._root in {a.value for a in self._onto.graph.get_ancestors(tid.identifier)}

    def closure(self, curies: Iterable[str]) -> Set[str]:
        """Ancestor closure of a whole annotation set."""
        out: Set[str] = set()
        for c in curies:
            out |= self.ancestors(c, include_self=True)
        return out


# --------------------------------------------------------------------------
# the IC model
# --------------------------------------------------------------------------
@dataclass
class ICModel:
    """A frozen frequency table and the IC values derived from it."""

    ic: Dict[str, float]
    freq: Dict[str, int]
    n_entities: int
    corpus: str
    normalized: bool
    smoothing: float
    max_ic: float
    unscoreable: Set[str] = field(default_factory=set)

    def __getitem__(self, curie: str) -> float:
        return self.ic[curie]

    def get(self, curie: str, default=None):
        return self.ic.get(curie, default)

    def describe(self) -> str:
        return (
            f"ICModel(corpus={self.corpus!r}, n_entities={self.n_entities}, "
            f"terms={len(self.ic)}, normalized={self.normalized}, "
            f"smoothing={self.smoothing}, max_ic={self.max_ic:.4f})"
        )


def build_ic_model(
    entity_annotations: Sequence[Iterable[str]],
    index: AncestorIndex,
    *,
    corpus: str = "pooled",
    normalize: bool = True,
    smoothing: float = 0.0,
) -> ICModel:
    """Build a Resnik IC table from a corpus of annotated entities.

    Parameters
    ----------
    entity_annotations
        One iterable of directly-annotated HPO CURIEs per corpus entity
        (a patient, or a disease). Order is irrelevant; duplicates within an
        entity are collapsed, because frequency is counted per *entity*, not
        per annotation.
    index
        AncestorIndex over the fixed scoring graph.
    corpus
        Label recorded on the model, one of CORPUS_CHOICES. Purely descriptive
        here -- the caller assembles the entity list.
    normalize
        Divide IC by ln(|C|), the IC of a term seen exactly once, putting values
        on [0, 1] with the root at 0. **Strongly recommended** whenever cohorts
        of different size are compared, since an unnormalised Resnik IC scales
        with ln(|C|) and would otherwise make the larger cohort look more
        informative for purely arithmetic reasons.
    smoothing
        Add-`smoothing` (Laplace) count applied to every term, so that terms
        absent from the corpus receive a finite IC instead of infinity. Only
        needed for `hpoa_replaced`, where a scored term may be missing from the
        background. Leave at 0 when the corpus is built from the union of the
        two arms, since every scored term then appears at least once.

    Notes
    -----
    **One frozen table for both arms.** Build the model once, from a corpus that
    covers both the aged and the curated annotations, then score both arms with
    it. Building a separate table per arm would put the two arms in different
    probability spaces and the difference in mean IC would no longer be
    interpretable -- a term's IC would change simply because the corpus changed
    underneath it. This mirrors the fixed-graph argument from the intrinsic
    version, and is the single most important thing to keep right.
    """
    if corpus not in CORPUS_CHOICES:
        warnings.warn(f"corpus={corpus!r} is not one of {CORPUS_CHOICES}")

    n = len(entity_annotations)
    if n == 0:
        raise ValueError("empty corpus")

    freq: Dict[str, int] = {}
    unscoreable: Set[str] = set()

    for terms in entity_annotations:
        expanded: Set[str] = set()
        for c in terms:
            anc = index.ancestors(c, include_self=True)
            if not anc:
                unscoreable.add(c)
                continue
            expanded |= anc
        for t in expanded:
            freq[t] = freq.get(t, 0) + 1

    denom = n + smoothing
    ic: Dict[str, float] = {}
    for t, f in freq.items():
        p = (f + smoothing) / denom
        p = min(p, 1.0)
        ic[t] = -math.log(p)

    max_ic = -math.log((1.0 + smoothing) / denom)
    if normalize and max_ic > 0:
        ic = {t: v / max_ic for t, v in ic.items()}

    return ICModel(
        ic=ic,
        freq=freq,
        n_entities=n,
        corpus=corpus,
        normalized=normalize,
        smoothing=smoothing,
        max_ic=1.0 if normalize else max_ic,
        unscoreable=unscoreable,
    )


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def score_annotation_set(
    terms: Iterable[str],
    model: ICModel,
    index: AncestorIndex,
    *,
    missing: str = "drop",
) -> Tuple[float, float, int, List[str]]:
    """Score one patient's directly-annotated terms.

    Only the *direct* annotations are scored, not their closure -- summing over
    the closure would count every ancestor of every term and reward deep
    annotations twice.

    Returns
    -------
    (sum_ic, mean_ic, n_scored, dropped)
    """
    vals: List[float] = []
    dropped: List[str] = []
    for c in terms:
        pid = index.primary_id(c)
        if pid is None or pid not in model.ic:
            if missing == "raise":
                raise KeyError(f"{c} has no IC in this model")
            dropped.append(c)
            continue
        vals.append(model.ic[pid])

    if not vals:
        return 0.0, float("nan"), 0, dropped
    return float(sum(vals)), float(sum(vals) / len(vals)), len(vals), dropped


# --------------------------------------------------------------------------
# HPOA handling (for corpus="hpoa_replaced")
# --------------------------------------------------------------------------
def load_hpoa_disease_annotations(
    hpoa_path,
    *,
    aspect: str = "P",
    drop_negated: bool = True,
) -> Dict[str, Set[str]]:
    """Parse `phenotype.hpoa` into {disease_id: {hpo_curie, ...}}.

    Only the phenotypic-abnormality aspect is kept by default. NOT-qualified
    (negated) rows are dropped, matching how the observed arm is scored.
    """
    import csv

    out: Dict[str, Set[str]] = {}
    with open(hpoa_path, "r", encoding="utf-8") as fh:
        rows = (ln for ln in fh if not ln.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        for row in reader:
            if aspect and row.get("aspect") != aspect:
                continue
            if drop_negated and (row.get("qualifier") or "").strip().upper() == "NOT":
                continue
            db = row.get("database_id")
            hp = row.get("hpo_id")
            if not db or not hp:
                continue
            out.setdefault(db, set()).add(hp)
    return out


def build_hpoa_replaced_corpus(
    hpoa_annotations: Mapping[str, Set[str]],
    cohort_patient_terms: Sequence[Iterable[str]],
    target_disease_id: str,
    *,
    patients_as_entities: bool = True,
) -> List[Set[str]]:
    """Delete the target disease from the HPOA corpus and insert the cohort.

    `patients_as_entities=True` adds each patient as its own corpus entity,
    which preserves within-cohort variation in how often a term is used.
    Setting it False instead inserts a single pseudo-disease carrying the union
    of the cohort's terms, which is closer to what an HPOA row actually is.
    """
    entities: List[Set[str]] = [
        set(v) for k, v in hpoa_annotations.items() if k != target_disease_id
    ]
    if patients_as_entities:
        entities.extend(set(t) for t in cohort_patient_terms)
    else:
        union: Set[str] = set()
        for t in cohort_patient_terms:
            union |= set(t)
        if union:
            entities.append(union)
    return entities
