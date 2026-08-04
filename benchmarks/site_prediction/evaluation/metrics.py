"""Metrics for ligand-binding-site prediction benchmarks.

The primary distance here is DCA: predicted center to closest ligand atom.
DCC conventionally means center-to-center and must not label this calculation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from benchmarks.site_prediction.schema import BindingSiteLabel, Coordinate, PredictedSite


@dataclass(frozen=True)
class SiteHit:
    """Evaluation of one true ligand site against a top-k prediction prefix."""

    target_id: str
    true_site_id: str
    top_k: int
    hit: bool
    dca_a: float
    matched_site_id: str
    matched_rank: int | None

    @property
    def dcc_a(self) -> float:
        """Deprecated compatibility alias; this metric is DCA, not DCC."""

        return self.dca_a


def euclidean_distance(a: Coordinate, b: Coordinate) -> float:
    return math.sqrt(sum((a[index] - b[index]) ** 2 for index in range(3)))


def dca_to_ligand_atoms(center: Coordinate, label: BindingSiteLabel) -> float:
    """Distance from predicted center to the closest atom of one reference ligand."""

    return min(euclidean_distance(center, atom) for atom in label.ligand_atoms)


def dcc_to_ligand_atoms(center: Coordinate, label: BindingSiteLabel) -> float:
    """Deprecated compatibility alias for :func:`dca_to_ligand_atoms`."""

    return dca_to_ligand_atoms(center, label)


def sorted_predictions(predictions: Iterable[PredictedSite]) -> list[PredictedSite]:
    return sorted(predictions, key=lambda site: (site.rank, site.site_id))


def evaluate_site_topk(
    *,
    target_id: str,
    label: BindingSiteLabel,
    predictions: Iterable[PredictedSite],
    top_k: int,
    threshold_a: float = 4.0,
) -> SiteHit:
    """Evaluate one reference ligand site using the first ``top_k`` predictions."""

    considered = sorted_predictions(predictions)[:top_k]
    if not considered:
        return SiteHit(
            target_id=target_id,
            true_site_id=label.site_id,
            top_k=top_k,
            hit=False,
            dca_a=float("inf"),
            matched_site_id="",
            matched_rank=None,
        )

    best_site = min(considered, key=lambda site: dca_to_ligand_atoms(site.center, label))
    best_dca = dca_to_ligand_atoms(best_site.center, label)
    return SiteHit(
        target_id=target_id,
        true_site_id=label.site_id,
        top_k=top_k,
        hit=best_dca <= threshold_a,
        dca_a=best_dca,
        matched_site_id=best_site.site_id,
        matched_rank=best_site.rank,
    )


def evaluate_target_standard(
    *,
    target_id: str,
    labels: Iterable[BindingSiteLabel],
    predictions: Iterable[PredictedSite],
    threshold_a: float = 4.0,
) -> dict[str, list[SiteHit]]:
    """
    Evaluate a target with the COACH420/HOLO4K-style Top-n and Top-(n+2) protocol.

    ``n`` is the number of true ligand sites in the target. DCA is measured from a
    predicted pocket center to the closest atom of the reference ligand.
    ``all_sites`` ignores ranking and asks whether the method's emitted portfolio
    contains the reference site at all. Correct ligand grouping is essential:
    the number of labels defines both ``n`` and the site-level denominator.
    """

    label_list = list(labels)
    n_sites = len(label_list)
    top_n = max(1, n_sites)
    top_n_plus_2 = top_n + 2
    prediction_list = sorted_predictions(predictions)
    all_sites = len(prediction_list)
    return {
        "top_n": [
            evaluate_site_topk(
                target_id=target_id,
                label=label,
                predictions=prediction_list,
                top_k=top_n,
                threshold_a=threshold_a,
            )
            for label in label_list
        ],
        "top_n_plus_2": [
            evaluate_site_topk(
                target_id=target_id,
                label=label,
                predictions=prediction_list,
                top_k=top_n_plus_2,
                threshold_a=threshold_a,
            )
            for label in label_list
        ],
        "all_sites": [
            evaluate_site_topk(
                target_id=target_id,
                label=label,
                predictions=prediction_list,
                top_k=all_sites,
                threshold_a=threshold_a,
            )
            for label in label_list
        ],
    }


def success_rate(hits: Iterable[SiteHit]) -> float:
    hit_list = list(hits)
    if not hit_list:
        return float("nan")
    return sum(1 for hit in hit_list if hit.hit) / len(hit_list)
