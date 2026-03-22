from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable

from ..events import TransformEvent
from .model import ActiveSiteEvent, ActiveSiteStateFrame, SiteReactionCouplingRow


def _flatten_strs(values: Iterable[Iterable[str] | str]) -> tuple[str, ...]:
    out: list[str] = []
    for value in values:
        if isinstance(value, str):
            out.append(value)
        else:
            out.extend(str(x) for x in value)
    return tuple(out)


def _event_atoms(evidence_bonds: list[dict]) -> set[int]:
    atoms: set[int] = set()
    for item in evidence_bonds:
        try:
            i = item.get("i", None)
            j = item.get("j", None)
        except Exception:
            continue
        if i is not None:
            atoms.add(int(i))
        if j is not None:
            atoms.add(int(j))
    return atoms


def build_site_reaction_couplings(
    transform_events: list[TransformEvent],
    site_states: list[ActiveSiteStateFrame],
    site_events: list[ActiveSiteEvent],
) -> list[SiteReactionCouplingRow]:
    states_by_frame_site: dict[tuple[int, str], ActiveSiteStateFrame] = {
        (row.frame, row.site_id): row for row in site_states
    }
    site_ids_by_frame: dict[int, list[str]] = defaultdict(list)
    for row in site_states:
        site_ids_by_frame[int(row.frame)].append(row.site_id)
    site_event_map: dict[tuple[int, str], ActiveSiteEvent] = {
        (row.frame, row.site_id): row for row in site_events
    }

    transforms_by_frame: dict[int, list[TransformEvent]] = defaultdict(list)
    for te in transform_events:
        transforms_by_frame[int(te.frame)].append(te)

    frames = sorted(
        set(transforms_by_frame)
        | {frame for frame, _site_id in states_by_frame_site}
        | {frame for frame, _site_id in site_event_map}
    )

    rows: list[SiteReactionCouplingRow] = []

    for frame in frames:
        curr_sites = sorted(site_ids_by_frame.get(frame, []))
        for site_id in curr_sites:
            curr = states_by_frame_site[(frame, site_id)]
            prev = states_by_frame_site.get((frame - 1, site_id))
            se = site_event_map.get((frame, site_id))

            site_changed = bool(se is not None)
            site_atoms = set(curr.core_members) | set(curr.state_members)
            if prev is not None:
                site_atoms |= set(prev.core_members) | set(prev.state_members)

            prev_assoc = tuple() if prev is None else tuple(prev.associated_species_labels)
            curr_assoc = tuple(curr.associated_species_labels)

            matched_transforms: list[TransformEvent] = []
            strengths: list[str] = []
            evidence_rows: list[dict] = []

            for te in transforms_by_frame.get(frame, []):
                te_atoms = _event_atoms(te.evidence_bonds)
                react_prod = set(te.reactants) | set(te.products)
                assoc_overlap = react_prod.intersection(set(prev_assoc) | set(curr_assoc))
                atom_overlap = te_atoms.intersection(site_atoms)

                if atom_overlap:
                    matched_transforms.append(te)
                    strengths.append("strong")
                    evidence_rows.extend(te.evidence_bonds)
                    continue

                if assoc_overlap:
                    matched_transforms.append(te)
                    strengths.append("medium")
                    evidence_rows.extend(te.evidence_bonds)
                    continue

            has_transform = bool(matched_transforms)
            if not has_transform and not site_changed:
                continue

            if has_transform and site_changed:
                coupling_type = "reaction_with_site_change"
            elif has_transform:
                coupling_type = "reaction_only_on_site"
            else:
                coupling_type = "site_change_only"

            link_strength = "none"
            if strengths:
                if "strong" in strengths:
                    link_strength = "strong"
                elif "medium" in strengths:
                    link_strength = "medium"

            rows.append(
                SiteReactionCouplingRow(
                    frame=frame,
                    site_id=site_id,
                    site_family=curr.site_family,
                    site_state_before="" if prev is None else prev.state_label,
                    site_state_after=curr.state_label,
                    site_changed=site_changed,
                    has_transform=has_transform,
                    transform_event_ids=tuple(te.event_id for te in matched_transforms),
                    transform_types=tuple(te.type for te in matched_transforms),
                    transform_reactants=_flatten_strs(te.reactants for te in matched_transforms),
                    transform_products=_flatten_strs(te.products for te in matched_transforms),
                    associated_species_before=prev_assoc if se is None else tuple(se.associated_species_before),
                    associated_species_after=curr_assoc if se is None else tuple(se.associated_species_after),
                    coupling_type=coupling_type,
                    link_strength=link_strength,
                    evidence_bonds=evidence_rows if evidence_rows else ([] if se is None else list(se.evidence_bonds)),
                )
            )

    return rows
