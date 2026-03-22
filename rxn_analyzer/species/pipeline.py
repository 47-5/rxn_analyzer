from __future__ import annotations
from collections import Counter
import numpy as np
from ase import Atoms

from .chem import (
    connected_components,
    wl_hash,
    formula,
    normalize_fragment_smiles,
    surface_signature,
    is_suspicious_smiles,
    SmilesStrategy,
)
from ..sites import SiteDefinition, SiteAssignment
from .model import SpeciesFrameSnapshot


class SpeciesLabeler:
    def __init__(
        self,
        wl_iters: int,
        ads_signature_mode: str,
        smiles_recompute_mode: str,
        smiles_strategy: SmilesStrategy,
        smiles_fallback_to_formula_if_suspicious: bool,
        smiles_allow_charged: bool,
        smiles_allow_dot: bool,
        site_signature_mode: str = "none",             # NEW
        site_definition: SiteDefinition | None = None, # NEW
    ):
        self.wl_iters = wl_iters
        self.ads_signature_mode = ads_signature_mode
        self.smiles_recompute_mode = smiles_recompute_mode
        self.smiles_strategy = smiles_strategy
        self.fallback_if_suspicious = smiles_fallback_to_formula_if_suspicious
        self.allow_charged_smiles = smiles_allow_charged
        self.allow_dot_smiles = smiles_allow_dot
        self.site_signature_mode = site_signature_mode
        self.site_definition = site_definition
        self._compkey_to_smiles: dict[tuple[frozenset[int], str], str] = {}

    def _compute_smiles(self, atoms: Atoms, comp: list[int], cov_edges: set[tuple[int, int]]) -> str:
        smi = self.smiles_strategy.compute(atoms, comp, cov_edges) or ""
        smi = normalize_fragment_smiles(atoms, comp, cov_edges, smi)
        if self.fallback_if_suspicious and smi and is_suspicious_smiles(
            smi,
            allow_charged=self.allow_charged_smiles,
            allow_dot=self.allow_dot_smiles,
        ):
            return ""
        return smi

    def _base_label(self, atoms: Atoms, comp: list[int], cov_edges: set[tuple[int, int]]) -> tuple[str, str]:
        f = formula(atoms, comp)
        wl = wl_hash(atoms, comp, cov_edges, iters=self.wl_iters)

        if self.smiles_recompute_mode == "on_change":
            key = (frozenset(comp), wl)
            smi = self._compkey_to_smiles.get(key)
            if smi is None:
                smi = self._compute_smiles(atoms, comp, cov_edges)
                self._compkey_to_smiles[key] = smi
        elif self.smiles_recompute_mode == "always":
            smi = self._compute_smiles(atoms, comp, cov_edges)
        else:
            raise ValueError(f"Unknown smiles_recompute_mode: {self.smiles_recompute_mode}")

        if smi:
            return f"{f}|smiles={smi}", wl
        return f"{f}|wl={wl}", wl

    def _format_site_suffix(self, asn: SiteAssignment) -> str:
        mode = (self.site_signature_mode or "none").strip().lower()
        if mode == "none":
            return ""

        if mode == "type":
            if asn.ambiguous:
                types = sorted({t for _sid, t in asn.candidate_sites})
                return f"|site_type={'+'.join(types)}|site_amb=1"
            return f"|site_type={asn.primary_site_type}"

        if mode == "id":
            if asn.ambiguous:
                ids = [sid for sid, _t in asn.candidate_sites]
                return f"|site={'+'.join(ids)}|site_amb=1"
            return f"|site={asn.primary_site_id}"

        if mode in ("type+id", "detailed"):
            if asn.ambiguous:
                pairs = [f"{sid}:{t}" for sid, t in asn.candidate_sites]
                return f"|site={'+'.join(pairs)}|site_amb=1"
            return f"|site={asn.primary_site_id}:{asn.primary_site_type}"

        raise ValueError(f"Unknown site_signature_mode: {self.site_signature_mode}")

    def _apply_site_signature(
        self,
        base: str,
        atoms: Atoms,
        comp: list[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
        slab_edges: set[tuple[int, int]],
    ) -> str:
        if self.site_definition is None:
            return base
        if (self.site_signature_mode or "none").strip().lower() == "none":
            return base

        asn = self.site_definition.assign_component(
            comp,
            slab_mask,
            ads_edges,
            atoms=atoms,
            slab_edges=slab_edges,
        )
        if asn is None:
            return base

        return base + self._format_site_suffix(asn)

    def _apply_ads_signature(
        self,
        base: str,
        atoms: Atoms,
        comp: list[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
        slab_edges: set[tuple[int, int]],
    ) -> str:
        surf = surface_signature(atoms, comp, slab_mask, ads_edges)

        if self.ads_signature_mode == "none":
            out = base
        elif self.ads_signature_mode == "coarse":
            out = f"{base}|ads" if surf is not None else base
        elif self.ads_signature_mode == "detailed":
            if surf is None:
                out = base
            else:
                out = (
                    f"{base}|ads(n={surf.n_ads_bonds},coord={list(surf.slab_coord_hist)},"
                    f"slab={list(surf.slab_elements)})"
                )
        elif self.ads_signature_mode == "topology":
            if surf is None:
                out = base
            else:
                pair_text = ",".join(surf.ads_pair_labels)
                out = f"{base}|ads(n={surf.n_ads_bonds},pairs=[{pair_text}])"
        else:
            raise ValueError(f"Unknown ads_signature_mode: {self.ads_signature_mode}")

        # NEW: append site info after ads signature
        out = self._apply_site_signature(out, atoms, comp, slab_mask, ads_edges, slab_edges)
        return out

    def labels_for_components(
        self,
        atoms: Atoms,
        components: list[list[int]],
        cov_edges: set[tuple[int, int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
        slab_edges: set[tuple[int, int]],
    ) -> list[str]:
        labels: list[str] = []
        for comp in components:
            base, _wl = self._base_label(atoms, comp, cov_edges)
            lbl = self._apply_ads_signature(base, atoms, comp, slab_mask, ads_edges, slab_edges)
            labels.append(lbl)
        return labels

    def ads_pairs_for_components(
        self,
        atoms: Atoms,
        components: list[list[int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
    ) -> dict[frozenset[int], list[str]]:
        out: dict[frozenset[int], list[str]] = {}
        for comp in components:
            surf = surface_signature(atoms, comp, slab_mask, ads_edges)
            out[frozenset(comp)] = list(surf.ads_pair_labels) if surf is not None else []
        return out

    def site_assignments_for_components(
        self,
        atoms: Atoms,
        components: list[list[int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
        slab_edges: set[tuple[int, int]],
    ) -> dict[frozenset[int], dict | None]:
        out: dict[frozenset[int], dict | None] = {}
        for comp in components:
            comp_key = frozenset(comp)
            if self.site_definition is None:
                out[comp_key] = None
                continue

            asn = self.site_definition.assign_component(
                comp,
                slab_mask,
                ads_edges,
                atoms=atoms,
                slab_edges=slab_edges,
            )
            if asn is None:
                out[comp_key] = None
                continue

            out[comp_key] = {
                "primary_site_id": asn.primary_site_id,
                "primary_site_type": asn.primary_site_type,
                "ambiguous": bool(asn.ambiguous),
                "candidate_sites": [list(x) for x in asn.candidate_sites],
                "n_ads_bonds": int(asn.n_ads_bonds),
                "touched_host_atoms": list(asn.touched_host_atoms),
            }
        return out

class SpeciesPipeline:
    def __init__(self, labeler: SpeciesLabeler):
        self.labeler = labeler

    def analyze(
        self,
        atoms: Atoms,
        cov_edges: set[tuple[int, int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
        slab_edges: set[tuple[int, int]],
    ) -> SpeciesFrameSnapshot:
        non_slab = {i for i in range(len(atoms)) if not slab_mask[i]}
        components = connected_components(len(atoms), cov_edges, non_slab)

        labels = self.labeler.labels_for_components(
            atoms,
            components,
            cov_edges,
            ads_edges,
            slab_mask,
            slab_edges,
        )
        multiset = Counter(labels)
        component_labels = {frozenset(component): label for component, label in zip(components, labels)}
        component_ads_pairs = self.labeler.ads_pairs_for_components(atoms, components, ads_edges, slab_mask)
        component_geometric_site_assignments = self.labeler.site_assignments_for_components(
            atoms,
            components,
            ads_edges,
            slab_mask,
            slab_edges,
        )

        return SpeciesFrameSnapshot(
            components=components,
            labels=labels,
            multiset=multiset,
            component_labels=component_labels,
            component_ads_pairs=component_ads_pairs,
            component_geometric_site_assignments=component_geometric_site_assignments,
        )
