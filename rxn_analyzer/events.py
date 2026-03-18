from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class TransformEvent:
    event_id: int
    frame: int
    type: str
    reactants: list[str]
    products: list[str]
    delta_counts: dict[str, int]
    evidence_bonds: list[dict] = field(default_factory=list)
    confidence: float | None = None