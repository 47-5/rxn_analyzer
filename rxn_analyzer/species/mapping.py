"""
把“上一帧的连通分量”与“当前帧的连通分量”做匹配/分块，并据此生成拆分计划（SplitPlan），用于描述反应物/产物组合与证据键
"""
from __future__ import annotations
from dataclasses import dataclass
from collections import Counter

from ..edges import EdgeType
from ..tracking import BondEvent


@dataclass
class SplitPlan:
    reactants: list[str]
    products: list[str]
    delta_counts: dict[str, int]
    evidence_bonds: list[dict]


class ComponentMapper:
    def __init__(
        self,
        *,
        split_by_component: bool,  # 是否启用
        component_overlap_min: float,  # 最小重叠比例
        component_overlap_mode: str,
        require_component_mapping: bool,  # 是否要求前后帧都必须有分量标签，否则直接返回 None
        keep_bond_evidence: bool,  # 是否分配 bond evidence
    ):
        self.split_by_component = split_by_component
        self.component_overlap_min = component_overlap_min
        self.component_overlap_mode = component_overlap_mode
        self.require_component_mapping = require_component_mapping
        self.keep_bond_evidence = keep_bond_evidence

        self._prev_comp_labels: dict[frozenset[int], str] | None = None

    def set_prev(self, comp_labels: dict[frozenset[int], str]) -> None:
        self._prev_comp_labels = comp_labels

    def _overlap_score(self, a: frozenset[int], b: frozenset[int]) -> float:
        inter = len(a & b)
        if inter == 0:
            return 0.0
        if self.component_overlap_mode == "old_to_new":
            return inter / max(1, len(a))
        if self.component_overlap_mode == "new_to_old":
            return inter / max(1, len(b))
        return min(inter / max(1, len(a)), inter / max(1, len(b)))

    def split(
        self,
        curr_comp_labels: dict[frozenset[int], str],
        bond_events: list[BondEvent],
    ) -> list[SplitPlan] | None:
        if not self.split_by_component:
            return None
        if self._prev_comp_labels is None:
            return None
        if self.require_component_mapping and (not self._prev_comp_labels or not curr_comp_labels):
            return None

        # --- 构建前后分量的“重叠邻接关系” ---
        prev_ids = list(self._prev_comp_labels.keys())
        curr_ids = list(curr_comp_labels.keys())

        adj_prev: dict[frozenset[int], set[frozenset[int]]] = {p: set() for p in prev_ids}  # adj_prev[p] = 能与旧分量 p 匹配的新分量集合
        adj_curr: dict[frozenset[int], set[frozenset[int]]] = {c: set() for c in curr_ids}  # adj_curr[c] = 能与新分量 c 匹配的旧分量集合

        for p in prev_ids:
            for c in curr_ids:
                sc = self._overlap_score(p, c)
                if sc > 0.0 and sc >= float(self.component_overlap_min):  # 仅当重叠比例 ≥ component_overlap_min
                    adj_prev[p].add(c)
                    adj_curr[c].add(p)

        has_any = any(adj_prev[p] for p in prev_ids)
        if not has_any:  # 如果没有任何匹配则返回
            return None

        # --- 找“连通块”（mapping block） ---
        visited_prev: set[frozenset[int]] = set()
        visited_curr: set[frozenset[int]] = set()
        blocks: list[tuple[set[frozenset[int]], set[frozenset[int]]]] = []

        for p0 in prev_ids:
            if p0 in visited_prev:
                continue
            if not adj_prev[p0]:  # 若某个旧分量没有任何新分量相连，则形成块：({p0}, set())
                visited_prev.add(p0)
                blocks.append(({p0}, set()))
                continue

            stack_prev = [p0]  # 否则，从 p0 出发穿越“旧↔新”图，收集所有连通的旧分量与新分量
            block_prev: set[frozenset[int]] = set()
            block_curr: set[frozenset[int]] = set()

            while stack_prev:
                p = stack_prev.pop()
                if p in visited_prev:
                    continue
                visited_prev.add(p)
                block_prev.add(p)
                for c in adj_prev[p]:
                    if c not in visited_curr:
                        visited_curr.add(c)
                        block_curr.add(c)
                        for p2 in adj_curr.get(c, set()):
                            if p2 not in visited_prev:
                                stack_prev.append(p2)

            blocks.append((block_prev, block_curr))

        for c0 in curr_ids:  # 之后把没有被访问到的新分量补上
            if c0 in visited_curr:
                continue
            visited_curr.add(c0)
            blocks.append((set(), {c0}))

        # --- 将每个 block 转成拆分计划 ---
        plans: list[SplitPlan] = []

        for prev_set, curr_set in blocks:
            reactants: list[str] = []
            products: list[str] = []
            for p in prev_set:
                reactants.append(self._prev_comp_labels[p])
            for c in curr_set:
                products.append(curr_comp_labels[c])

            if not reactants and not products:
                continue
            if not reactants or not products:
                continue

            if Counter(reactants) == Counter(products):
                continue

            ev_bonds: list[dict] = []
            if self.keep_bond_evidence:
                prev_atoms = set().union(*prev_set) if prev_set else set()
                curr_atoms = set().union(*curr_set) if curr_set else set()
                atoms_union = prev_atoms | curr_atoms

                for be in bond_events:
                    if be.edge.type == EdgeType.SLAB:
                        continue
                    if be.edge.type == EdgeType.COV:
                        if (be.edge.i in atoms_union) and (be.edge.j in atoms_union):
                            ev_bonds.append(
                                {
                                    "edge_type": be.edge.type.value,
                                    "action": be.action,
                                    "i": be.edge.i,
                                    "j": be.edge.j,
                                    "distance": be.evidence.distance,
                                    "threshold_form": be.evidence.threshold_form,
                                    "threshold_break": be.evidence.threshold_break,
                                }
                            )
                    elif be.edge.type == EdgeType.ADS:
                        if (be.edge.i in atoms_union) or (be.edge.j in atoms_union):
                            ev_bonds.append(
                                {
                                    "edge_type": be.edge.type.value,
                                    "action": be.action,
                                    "i": be.edge.i,
                                    "j": be.edge.j,
                                    "distance": be.evidence.distance,
                                    "threshold_form": be.evidence.threshold_form,
                                    "threshold_break": be.evidence.threshold_break,
                                }
                            )

            delta_counts: dict[str, int] = {}
            for r in reactants:
                delta_counts[r] = delta_counts.get(r, 0) - 1
            for p in products:
                delta_counts[p] = delta_counts.get(p, 0) + 1

            plans.append(
                SplitPlan(
                    reactants=reactants,
                    products=products,
                    delta_counts=delta_counts,
                    evidence_bonds=ev_bonds,
                )
            )

        return plans
