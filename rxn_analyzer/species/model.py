from __future__ import annotations

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class SpeciesFrameSnapshot:
    components: list[list[int]]
    labels: list[str]
    multiset: Counter[str]
    component_labels: dict[frozenset[int], str]
    component_ads_pairs: dict[frozenset[int], list[str]]
    component_geometric_site_assignments: dict[frozenset[int], dict | None]
