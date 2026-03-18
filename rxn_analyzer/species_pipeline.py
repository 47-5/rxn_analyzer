from __future__ import annotations
from dataclasses import dataclass
from collections import Counter
import numpy as np
from ase import Atoms

from .species import (
    connected_components,
    wl_hash,
    formula,
    surface_signature,
    is_suspicious_smiles,
    SmilesStrategy,
)
from .sites import SiteDefinition, SiteAssignment


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
        comp: list[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
    ) -> str:
        if self.site_definition is None:
            return base
        if (self.site_signature_mode or "none").strip().lower() == "none":
            return base

        asn = self.site_definition.assign_component(comp, slab_mask, ads_edges)
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
                    f"{base}|ads(n={surf.n_ads_bonds},coord={list(surf.slab_coord_hist)},slab={list(surf.slab_elements)})"
                )
        else:
            raise ValueError(f"Unknown ads_signature_mode: {self.ads_signature_mode}")

        # NEW: append site info after ads signature
        out = self._apply_site_signature(out, comp, slab_mask, ads_edges)
        return out

    def labels_for_components(
        self,
        atoms: Atoms,
        components: list[list[int]],
        cov_edges: set[tuple[int, int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
    ) -> list[str]:
        labels: list[str] = []
        for comp in components:
            base, _wl = self._base_label(atoms, comp, cov_edges)
            lbl = self._apply_ads_signature(base, atoms, comp, slab_mask, ads_edges)
            labels.append(lbl)
        return labels


@dataclass
class SpeciesFrame:
    comps: list[list[int]]
    labels: list[str]
    multiset: Counter[str]
    comp_labels: dict[frozenset[int], str]


class SpeciesPipeline:
    def __init__(self, labeler: SpeciesLabeler):
        self.labeler = labeler

    def analyze(
        self,
        atoms: Atoms,
        cov_edges: set[tuple[int, int]],
        ads_edges: set[tuple[int, int]],
        slab_mask: np.ndarray,
    ) -> SpeciesFrame:
        non_slab = {i for i in range(len(atoms)) if not slab_mask[i]}
        comps = connected_components(len(atoms), cov_edges, non_slab)

        labels = self.labeler.labels_for_components(atoms, comps, cov_edges, ads_edges, slab_mask)
        multiset = Counter(labels)
        comp_labels = {frozenset(c): lbl for c, lbl in zip(comps, labels)}

        return SpeciesFrame(
            comps=comps,
            labels=labels,
            multiset=multiset,
            comp_labels=comp_labels,
        )