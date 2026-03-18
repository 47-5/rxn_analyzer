from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, Set, List

from .edges import Edge, EdgeEvidence, EdgeType


@dataclass
class EdgeTrack:
    """
    State machine for a single edge.
    - confirmed: current confirmed state (0/1) 当前确认状态（0/1），才是系统“承认”的键状态。
    - candidate: proposed new state 当前正在尝试切换到的状态（0/1），None 表示没有切换趋势。
    - streak: consecutive frames candidate observed          candidate 连续观察到的帧数。
    - cooldown: frames remaining before allowing another switch  切换后冷却期（阻止快速翻转）。
    """
    confirmed: int = 0
    candidate: int | None = None
    streak: int = 0
    cooldown: int = 0


@dataclass(frozen=True, slots=True)
class BondEvent:
    event_id: int
    frame: int
    edge: Edge
    action: str          # "form" / "break"
    evidence: EdgeEvidence


class EdgeTracker:
    """
    Converts instantaneous (hysteresis-processed) edge states into confirmed states
    using:
      - persist: number of consecutive frames needed to switch
      - cooldown: lockout after switching to prevent oscillations
    """
    def __init__(self, persist: int = 5, cooldown: int = 3):
        if persist < 1:
            raise ValueError("persist must be >= 1")
        if cooldown < 0:
            raise ValueError("cooldown must be >= 0")
        self.persist = persist
        self.cooldown = cooldown
        self.tracks: Dict[Edge, EdgeTrack] = {}  # Edge → EdgeTrack 的表。

    def get_confirmed(self, e: Edge) -> int:
        # 获取某条边“确认状态”，若不存在就默认 0。
        t = self.tracks.get(e)
        return 0 if t is None else t.confirmed  # 给 build_inst_state_strict_hysteresis 用。不存在的边默认为 confirmed = 0。

    def update(
        self,
        frame: int,
        inst_state: Dict[Edge, Tuple[int, EdgeEvidence]],  # 当前帧“瞬时状态” (严格迟滞判定过的 0/1)
        universe: Set[Edge],
        next_event_id: int,
    ) -> tuple[List[BondEvent], int]:
        """
        inst_state must already reflect STRICT hysteresis using previous confirmed state.
        For edges in universe but absent in inst_state, treat instantaneous=0.前提：inst_state 已经是严格迟滞后的瞬时状态。
        """
        events: List[BondEvent] = []

        # ensure tracks for all edges in universe  确保每条边都应有状态机，即便它之前从未出现过。
        for e in universe:
            if e not in self.tracks:
                self.tracks[e] = EdgeTrack(confirmed=0)

        for e in universe:
            x, ev = inst_state.get(e, (0, EdgeEvidence()))  # x 是瞬时状态（0/1）  若没有 inst_state 记录，使用 (0, EdgeEvidence())。
            t = self.tracks[e]

            if t.cooldown > 0:  # 冷却期内 不允许切换。直接清空 candidate/streak，防止累计
                t.cooldown -= 1
                # cooldown: block switching; also clear candidate memory
                t.candidate = None
                t.streak = 0
                continue

            if x == t.confirmed:  # 没有切换趋势（当前状态和已确认状态相同），候选清空
                t.candidate = None
                t.streak = 0
                continue

            # propose switching
            if t.candidate is None or t.candidate != x:  # 如果候选状态变化了，重新开始计数
                t.candidate = x
                t.streak = 1
            else:  # 如果保持一致，就 streak++
                t.streak += 1

            if t.streak >= self.persist:
                old = t.confirmed
                t.confirmed = x
                t.candidate = None
                t.streak = 0
                t.cooldown = self.cooldown

                action = "form" if (old == 0 and x == 1) else "break"
                events.append(BondEvent(next_event_id, frame, e, action, ev))
                next_event_id += 1

        return events, next_event_id

    def confirmed_edges(self, edge_type: EdgeType) -> Set[tuple[int, int]]:  # 按类型筛选确认 ON 的边。 返回的是 (i, j) 元组，保证 i≤j。

        out: Set[tuple[int, int]] = set()
        for e, t in self.tracks.items():
            if e.type == edge_type and t.confirmed == 1:
                out.add((e.i, e.j) if e.i <= e.j else (e.j, e.i))
        return out