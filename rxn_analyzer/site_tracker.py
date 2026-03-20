from __future__ import annotations

from collections import Counter

from .site_model import JointReactiveSiteReaction, ReactiveSiteEvent, ReactiveSiteStateFrame
from .tracking import BondEvent


def _bond_event_to_dict(be: BondEvent) -> dict[str, object]:
    return {
        "edge_type": be.edge.type.value,
        "action": be.action,
        "i": be.edge.i,
        "j": be.edge.j,
        "distance": be.evidence.distance,
        "threshold_form": be.evidence.threshold_form,
        "threshold_break": be.evidence.threshold_break,
    }


class ReactiveSiteTracker:
    def __init__(self) -> None:
        self.prev_states: dict[str, ReactiveSiteStateFrame] = {}
        self.events: list[ReactiveSiteEvent] = []
        self.joint_reactions: list[JointReactiveSiteReaction] = []

    @staticmethod
    def _filter_evidence(site_prev: ReactiveSiteStateFrame | None,
                         site_curr: ReactiveSiteStateFrame,
                         bond_events: list[BondEvent]) -> list[dict[str, object]]:
        atoms_of_interest = set(site_curr.core_members) | set(site_curr.state_members)
        if site_prev is not None:
            atoms_of_interest |= set(site_prev.core_members) | set(site_prev.state_members)

        out: list[dict[str, object]] = []
        for be in bond_events:
            if be.edge.i in atoms_of_interest or be.edge.j in atoms_of_interest:
                out.append(_bond_event_to_dict(be))
        return out

    def update(
        self,
        frame: int,
        states: list[ReactiveSiteStateFrame],
        bond_events: list[BondEvent],
        next_event_id: int,
    ) -> int:
        curr_map = {s.site_id: s for s in states}

        for site_id, curr in curr_map.items():
            prev = self.prev_states.get(site_id)
            if prev is None:
                continue

            prev_members = set(prev.state_members)
            curr_members = set(curr.state_members)
            added = tuple(sorted(curr_members - prev_members))
            removed = tuple(sorted(prev_members - curr_members))
            state_changed = prev.state_label != curr.state_label
            species_changed = Counter(prev.associated_species_labels) != Counter(curr.associated_species_labels)

            if not state_changed and not added and not removed and not species_changed:
                continue

            evidence = self._filter_evidence(prev, curr, bond_events)

            self.events.append(
                ReactiveSiteEvent(
                    event_id=next_event_id,
                    frame=frame,
                    site_id=site_id,
                    event_type="site_state_change",
                    old_state=prev.state_label,
                    new_state=curr.state_label,
                    added_members=added,
                    removed_members=removed,
                    bound_species_before=tuple(prev.associated_species_labels),
                    bound_species_after=tuple(curr.associated_species_labels),
                    evidence_bonds=evidence,
                )
            )
            next_event_id += 1

            if species_changed or state_changed:
                self.joint_reactions.append(
                    JointReactiveSiteReaction(
                        event_id=next_event_id,
                        frame=frame,
                        site_id=site_id,
                        site_state_before=prev.state_label,
                        site_state_after=curr.state_label,
                        reactants=tuple(prev.associated_species_labels),
                        products=tuple(curr.associated_species_labels),
                        bound_species_before=tuple(prev.associated_species_labels),
                        bound_species_after=tuple(curr.associated_species_labels),
                        evidence_bonds=evidence,
                    )
                )
                next_event_id += 1

        self.prev_states = curr_map
        return next_event_id
