from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from ase import Atoms

from .site_model import ReactiveSiteDefinition, ReactiveSiteStateFrame
from .site_rules import formula_from_atoms, topology_from_state_components


@dataclass(frozen=True)
class ReactiveSiteFrameBatch:
    states: list[ReactiveSiteStateFrame]


class ReactiveSitePipeline:
    def __init__(self, definition: ReactiveSiteDefinition):
        self.definition = definition

    @staticmethod
    def _component_lookup(components: list[list[int]]) -> tuple[dict[int, tuple[int, ...]], dict[tuple[int, ...], str]]:
        atom_to_component: dict[int, tuple[int, ...]] = {}
        component_to_label: dict[tuple[int, ...], str] = {}
        for comp in components:
            key = tuple(sorted(comp))
            for atom in comp:
                atom_to_component[atom] = key
        return atom_to_component, component_to_label

    @staticmethod
    def _build_cov_adjacency(cov_edges: set[tuple[int, int]]) -> dict[int, set[int]]:
        adj: dict[int, set[int]] = defaultdict(set)
        for i, j in cov_edges:
            adj[i].add(j)
            adj[j].add(i)
        return adj

    @staticmethod
    def _build_core_edge_map(
        site_core: frozenset[int],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
    ) -> dict[int, set[int]]:
        core_edge_map: dict[int, set[int]] = defaultdict(set)
        n = int(len(slab_mask))
        for i, j in ads_edges:
            if i in site_core and 0 <= j < n and not bool(slab_mask[j]):
                core_edge_map[j].add(i)
            elif j in site_core and 0 <= i < n and not bool(slab_mask[i]):
                core_edge_map[i].add(j)
        return core_edge_map

    @staticmethod
    def _component_keys_for_atoms(
        atoms: set[int],
        atom_to_component: dict[int, tuple[int, ...]],
    ) -> set[tuple[int, ...]]:
        out: set[tuple[int, ...]] = set()
        for atom in atoms:
            comp_key = atom_to_component.get(atom)
            if comp_key is not None:
                out.add(comp_key)
        return out

    @staticmethod
    def _build_core_contacts(
        site_core: frozenset[int],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
    ) -> set[int]:
        contacts: set[int] = set()
        n = int(len(slab_mask))
        for i, j in ads_edges:
            if i in site_core and 0 <= j < n and not bool(slab_mask[j]):
                contacts.add(j)
            elif j in site_core and 0 <= i < n and not bool(slab_mask[i]):
                contacts.add(i)
        return contacts

    def analyze(
        self,
        frame: int,
        atoms: Atoms,
        slab_mask: np.ndarray,
        cov_edges: set[tuple[int, int]],
        ads_edges: set[tuple[int, int]],
        components: list[list[int]],
        labels: list[str],
    ) -> ReactiveSiteFrameBatch:
        if not self.definition.sites:
            return ReactiveSiteFrameBatch(states=[])

        atom_to_component, component_to_label = self._component_lookup(components)
        for comp, label in zip(components, labels):
            component_to_label[tuple(sorted(comp))] = label

        cov_adj = self._build_cov_adjacency(cov_edges)
        states: list[ReactiveSiteStateFrame] = []

        n_atoms = len(atoms)
        syms = atoms.get_chemical_symbols()

        for site in self.definition.sites:
            invalid_core = [a for a in site.core_members if a < 0 or a >= n_atoms or not bool(slab_mask[a])]
            if invalid_core and self.definition.strict_core_validation:
                raise ValueError(
                    f"Reactive site '{site.site_id}' has core_members not in current slab mask: {sorted(invalid_core)}"
                )
            valid_core = frozenset(a for a in site.core_members if 0 <= a < n_atoms and bool(slab_mask[a]))

            core_edge_map = self._build_core_edge_map(valid_core, ads_edges, slab_mask)
            direct_contacts = self._build_core_contacts(valid_core, ads_edges, slab_mask)
            intrinsic_seed_members: set[int] = set()

            for atom in site.initial_state_members:
                if 0 <= atom < n_atoms and not bool(slab_mask[atom]):
                    intrinsic_seed_members.add(atom)

            intrinsic_component_keys = self._component_keys_for_atoms(intrinsic_seed_members, atom_to_component)
            associated_component_keys = self._component_keys_for_atoms(direct_contacts, atom_to_component)

            intrinsic_members: set[int] = set()
            for comp_key in intrinsic_component_keys:
                for atom in comp_key:
                    if bool(slab_mask[atom]):
                        continue
                    if site.allowed_state_elements and syms[atom] not in site.allowed_state_elements:
                        continue
                    intrinsic_members.add(atom)

            component_atom_cache: dict[tuple[int, ...], set[int]] = {}
            candidate_component_atoms: dict[tuple[int, ...], set[int]] = {}
            for comp_key in associated_component_keys:
                component_atoms: set[int] = set()
                for atom in comp_key:
                    if bool(slab_mask[atom]):
                        continue
                    if site.allowed_state_elements and syms[atom] not in site.allowed_state_elements:
                        continue
                    component_atoms.add(atom)
                if component_atoms:
                    component_atom_cache[comp_key] = component_atoms
                    candidate_component_atoms[comp_key] = component_atoms

            incorporated_members: set[int] = set()
            promoted_component_keys: set[tuple[int, ...]] = set()
            site_network_members = set(intrinsic_members)

            changed = True
            while changed:
                changed = False
                for comp_key, comp_atoms in sorted(candidate_component_atoms.items()):
                    if comp_key in promoted_component_keys:
                        continue

                    core_anchor_atoms: set[int] = set()
                    for atom in comp_atoms:
                        core_anchor_atoms.update(core_edge_map.get(atom, set()))

                    intrinsic_anchor_atoms: set[int] = set()
                    for atom in comp_atoms:
                        for nbr in cov_adj.get(atom, set()):
                            if nbr in site_network_members:
                                intrinsic_anchor_atoms.add(nbr)

                    all_anchor_atoms = set(core_anchor_atoms) | set(intrinsic_anchor_atoms)
                    score = 0
                    if core_anchor_atoms:
                        score += 1
                    if intrinsic_anchor_atoms:
                        score += 1
                    if len(all_anchor_atoms) >= 2:
                        score += 1

                    if score >= 2:
                        promoted_component_keys.add(comp_key)
                        incorporated_members.update(comp_atoms)
                        site_network_members.update(comp_atoms)
                        changed = True

            associated_members: set[int] = set()
            filtered_associated_component_keys = sorted(
                comp_key for comp_key in candidate_component_atoms if comp_key not in promoted_component_keys
            )
            for comp_key in filtered_associated_component_keys:
                associated_members.update(candidate_component_atoms[comp_key])

            state_members = intrinsic_members | incorporated_members | associated_members

            if site.max_state_members is not None and len(state_members) > site.max_state_members:
                seeds = set(direct_contacts) | set(intrinsic_seed_members)
                ordered = sorted(state_members, key=lambda a: (0 if a in seeds else 1, a))
                keep = set(ordered[: site.max_state_members])
                intrinsic_members &= keep
                incorporated_members &= keep
                associated_members &= keep
                state_members = keep

            filtered_component_keys: list[tuple[int, ...]] = []
            associated_species_labels: list[str] = []
            for comp_key in filtered_associated_component_keys:
                kept = tuple(sorted(a for a in comp_key if a in associated_members))
                if not kept:
                    continue
                filtered_component_keys.append(tuple(sorted(comp_key)))
                lbl = component_to_label.get(tuple(sorted(comp_key)))
                if lbl:
                    associated_species_labels.append(lbl)

            state_member_list = sorted(state_members)
            intrinsic_member_list = sorted(intrinsic_members)
            incorporated_member_list = sorted(incorporated_members)
            associated_member_list = sorted(associated_members)

            intrinsic_formula = formula_from_atoms(atoms, intrinsic_member_list)
            incorporated_formula = formula_from_atoms(atoms, incorporated_member_list)
            associated_formula = formula_from_atoms(atoms, associated_member_list)
            state_formula = formula_from_atoms(atoms, sorted(intrinsic_members | incorporated_members))

            def build_topo_components(member_list: list[int]) -> list[list[int]]:
                topo_components: list[list[int]] = []
                visited: set[int] = set()
                member_set = set(member_list)
                for atom in member_list:
                    if atom in visited:
                        continue
                    stack = [atom]
                    visited.add(atom)
                    comp = [atom]
                    while stack:
                        x = stack.pop()
                        for y in cov_adj.get(x, set()):
                            if y in member_set and y not in visited:
                                visited.add(y)
                                stack.append(y)
                                comp.append(y)
                    topo_components.append(sorted(comp))
                return topo_components

            intrinsic_topology = topology_from_state_components(
                atoms,
                build_topo_components(intrinsic_member_list),
                core_contacts=direct_contacts,
            )
            incorporated_topology = topology_from_state_components(
                atoms,
                build_topo_components(incorporated_member_list),
                core_contacts=direct_contacts,
            )
            associated_topology = topology_from_state_components(
                atoms,
                build_topo_components(associated_member_list),
                core_contacts=direct_contacts,
            )
            state_topology = topology_from_state_components(
                atoms,
                build_topo_components(sorted(intrinsic_members | incorporated_members)),
                core_contacts=direct_contacts,
            )

            base_label = state_topology if state_topology not in {"", "bare"} else state_formula
            associated_label = (
                associated_topology if associated_topology not in {"", "bare"} else associated_formula
            )
            if associated_label not in {"", "bare"}:
                state_label = f"{base_label if base_label not in {'', 'bare'} else 'bare'} + assoc({associated_label})"
            else:
                state_label = base_label if base_label not in {"", "bare"} else "bare"

            descriptors: dict[str, object] = {
                "n_core_members": len(valid_core),
                "n_intrinsic_members": len(intrinsic_member_list),
                "n_incorporated_members": len(incorporated_member_list),
                "n_associated_members": len(associated_member_list),
                "n_state_members": len(state_member_list),
                "n_associated_species": len(associated_species_labels),
                "direct_core_contacts": sorted(direct_contacts),
                "intrinsic_formula": intrinsic_formula,
                "intrinsic_topology": intrinsic_topology,
                "incorporated_formula": incorporated_formula,
                "incorporated_topology": incorporated_topology,
                "associated_formula": associated_formula,
                "associated_topology": associated_topology,
                "state_formula": state_formula,
                "state_topology": state_topology,
            }

            states.append(
                ReactiveSiteStateFrame(
                    frame=frame,
                    site_id=site.site_id,
                    site_family=site.site_family,
                    core_members=tuple(sorted(valid_core)),
                    intrinsic_members=tuple(intrinsic_member_list),
                    incorporated_members=tuple(incorporated_member_list),
                    associated_members=tuple(associated_member_list),
                    state_members=tuple(state_member_list),
                    intrinsic_formula=intrinsic_formula,
                    intrinsic_topology=intrinsic_topology,
                    incorporated_formula=incorporated_formula,
                    incorporated_topology=incorporated_topology,
                    associated_formula=associated_formula,
                    associated_topology=associated_topology,
                    state_formula=state_formula,
                    state_topology=state_topology,
                    state_label=state_label,
                    associated_species_labels=tuple(sorted(associated_species_labels)),
                    attached_components=tuple(filtered_component_keys),
                    descriptors=descriptors,
                )
            )

        return ReactiveSiteFrameBatch(states=states)
