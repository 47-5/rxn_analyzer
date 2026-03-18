from __future__ import annotations
from dataclasses import dataclass
from ase.data import covalent_radii


@dataclass(frozen=True)
class DistanceHysteresisParams:
    """
    Strict hysteresis:
      - if confirmed OFF: ON when d <= d_form
      - if confirmed ON : stay ON while d <= d_break  (break when d > d_break)
    """
    scale_form: float
    scale_break: float
    max_neighbors: int | None = None

    def thresholds(self, Zi: int, Zj: int) -> tuple[float, float]:
        ri = float(covalent_radii[Zi])
        rj = float(covalent_radii[Zj])
        base = ri + rj
        return self.scale_form * base, self.scale_break * base


@dataclass(frozen=True)
class Criteria:
    # 集中管理“所有边类型的成键参数”和“时间迟滞参数”。
    persist: int = 5  # 时间参数（给 EdgeTracker 用）
    cooldown: int = 3

    cov: DistanceHysteresisParams = DistanceHysteresisParams(1.15, 1.25, None)  # 距离参数（给成键判定用）
    ads: DistanceHysteresisParams = DistanceHysteresisParams(1.25, 1.40, None)
    slab: DistanceHysteresisParams = DistanceHysteresisParams(1.35, 1.55, 12)  # slab 的阈值更宽松 + 限制最大邻居数 因为金属表面原子常常配位很多，需要“剪枝”。