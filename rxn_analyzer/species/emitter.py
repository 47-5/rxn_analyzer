from __future__ import annotations

from ..events import TransformEvent
from ..edges import EdgeType
from ..tracking import BondEvent
from ..graph import add_transform_bipartite
from .mapping import SplitPlan


class EventIdCounter:
    def __init__(self, start: int = 1):
        self.value = int(start)

    def next_id(self) -> int:
        v = self.value
        self.value += 1
        return v

    def advance_to(self, new_value: int) -> None:
        self.value = int(new_value)


def classify_transform(delta: dict[str, int], bond_events: list[BondEvent]) -> str:
    changed_types = {be.edge.type for be in bond_events}
    has_species_change = any(v != 0 for v in delta.values())

    def _is_ads_label(s: str) -> bool:
        return "|ads" in s

    ads_delta = 0
    non_ads_delta = 0
    for sp, dv in delta.items():
        if dv == 0:
            continue
        if _is_ads_label(sp):
            ads_delta += dv
        else:
            non_ads_delta += dv

    has_ads_label_change = (ads_delta != 0)

    changed_types_noslab = {t for t in changed_types if t != EdgeType.SLAB}

    if not has_species_change and changed_types == {EdgeType.SLAB}:
        return "surface_reconstruct"

    if has_species_change:
        if has_ads_label_change and EdgeType.COV not in changed_types_noslab:
            return "adsorb/desorb/migrate"

        if EdgeType.COV not in changed_types_noslab and EdgeType.ADS in changed_types_noslab:
            return "adsorb/desorb/migrate"

        return "reaction"

    if EdgeType.COV in changed_types_noslab or EdgeType.ADS in changed_types_noslab:
        return "rearrange"
    return "surface_reconstruct" if EdgeType.SLAB in changed_types else "rearrange"


class TransformEmitter:
    def __init__(self, graph, id_counter: EventIdCounter):
        self.graph = graph
        self.id_counter = id_counter
        self.transform_events: list[TransformEvent] = []

    def emit_from_delta(
        self,
        frame: int,
        delta: dict[str, int],
        bond_events: list[BondEvent],
        *,
        ttype_override: str | None,
        add_to_graph: bool,  # 是否写入图结构
    ) -> None:
        reactants, products = [], []
        for k, dv in sorted(delta.items()):
            if dv < 0:
                reactants += [k] * (-dv)
            elif dv > 0:
                products += [k] * dv

        ev_bonds = [
            {
                "edge_type": be.edge.type.value,
                "action": be.action,
                "i": be.edge.i,
                "j": be.edge.j,
                "distance": be.evidence.distance,
                "threshold_form": be.evidence.threshold_form,
                "threshold_break": be.evidence.threshold_break,
            }
            for be in bond_events
        ]

        ttype = ttype_override or classify_transform(delta, bond_events)

        eid = self.id_counter.next_id()
        te = TransformEvent(
            event_id=eid,
            frame=frame,
            type=ttype,
            reactants=reactants,
            products=products,
            delta_counts=delta,
            evidence_bonds=ev_bonds,
            confidence=None,
        )
        self.transform_events.append(te)

        if add_to_graph:
            add_transform_bipartite(self.graph, eid, reactants, products, ttype)

    def emit_from_plans(
        self,
        frame: int,
        plans: list[SplitPlan],
        bond_events: list[BondEvent],
        *,
        add_to_graph: bool,
    ) -> int:
        emitted = 0
        for plan in plans:
            ttype = classify_transform(plan.delta_counts, bond_events)

            eid = self.id_counter.next_id()
            te = TransformEvent(
                event_id=eid,
                frame=frame,
                type=ttype,
                reactants=plan.reactants,
                products=plan.products,
                delta_counts=plan.delta_counts,
                evidence_bonds=plan.evidence_bonds,
                confidence=None,
            )
            self.transform_events.append(te)
            emitted += 1

            if add_to_graph:
                add_transform_bipartite(self.graph, eid, plan.reactants, plan.products, ttype)

        return emitted
