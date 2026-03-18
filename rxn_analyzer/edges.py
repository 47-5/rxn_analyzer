from __future__ import annotations
from dataclasses import dataclass
from enum import Enum


class EdgeType(str, Enum):
    COV = "covalent"       # non-slab ↔ non-slab
    ADS = "adsorption"     # slab ↔ non-slab
    SLAB = "slab"          # slab ↔ slab


@dataclass(frozen=True, slots=True)  # frozen=True 使对象不可变（hashable），可以作为 dict / set 的 key
class Edge:
    i: int
    j: int
    type: EdgeType

    def canonical(self) -> "Edge":
        return self if self.i <= self.j else Edge(self.j, self.i, self.type)  # 保证边的表示唯一化（i ≤ j）


@dataclass(frozen=True, slots=True)
class EdgeEvidence:
    # 记录“为什么判定成键 / 断键”的证据
    distance: float | None = None
    threshold_form: float | None = None
    threshold_break: float | None = None