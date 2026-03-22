from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ase import Atoms
from ase.io import read

from rxn_analyzer.config_loader import (
    RunOverrides,
    _build_criteria,
    _build_slab_def,
    _get_section,
    _resolve_path,
    apply_overrides,
    load_raw_config,
)
from rxn_analyzer.criteria import Criteria
from rxn_analyzer.edge_pipeline import build_inst_state_strict_hysteresis
from rxn_analyzer.edges import EdgeType
from rxn_analyzer.tracking import EdgeTracker


@dataclass(frozen=True)
class GeneratedActiveSite:
    site_id: str
    family: str
    core_members: tuple[int, ...]
    initial_state_members: tuple[int, ...]
    allowed_state_elements: tuple[str, ...]
    max_state_members: int | None
    rule_profile: str
    metadata: dict[str, Any]

    def to_mapping(self, *, output_index_base: int = 0) -> dict[str, Any]:
        shift = 0 if output_index_base == 0 else 1
        return {
            "id": self.site_id,
            "family": self.family,
            "core_members": [int(x) + shift for x in self.core_members],
            "initial_state_members": [int(x) + shift for x in self.initial_state_members],
            "allowed_state_elements": list(self.allowed_state_elements),
            "max_state_members": self.max_state_members,
            "rule_profile": self.rule_profile,
            "metadata": self.metadata,
        }


def _load_structure(path: str, frame_index: int) -> Atoms:
    return read(path, index=int(frame_index))


def _build_edges(
    atoms: Atoms,
    host_mask,
    criteria: Criteria,
) -> tuple[set[tuple[int, int]], set[tuple[int, int]], set[tuple[int, int]]]:
    tracker = EdgeTracker(persist=1, cooldown=0)
    inst, _universe = build_inst_state_strict_hysteresis(atoms, host_mask, criteria, tracker)
    cov_edges: set[tuple[int, int]] = set()
    ads_edges: set[tuple[int, int]] = set()
    host_edges: set[tuple[int, int]] = set()
    for edge, (on, _ev) in inst.items():
        if on != 1:
            continue
        ij = (edge.i, edge.j) if edge.i <= edge.j else (edge.j, edge.i)
        if edge.type == EdgeType.COV:
            cov_edges.add(ij)
        elif edge.type == EdgeType.ADS:
            ads_edges.add(ij)
        elif edge.type == EdgeType.SLAB:
            host_edges.add(ij)
    return cov_edges, ads_edges, host_edges


def _adjacency(edges: set[tuple[int, int]]) -> dict[int, set[int]]:
    adj: dict[int, set[int]] = defaultdict(set)
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _pick_best_framework_neighbor(atoms: Atoms, oxygen: int, candidates: list[int]) -> int:
    if len(candidates) == 1:
        return candidates[0]
    ox = atoms.positions[oxygen]
    return min(candidates, key=lambda idx: float(((atoms.positions[idx] - ox) ** 2).sum()))


def generate_bas_oh_bridge_sites(
    atoms: Atoms,
    host_mask,
    cov_edges: set[tuple[int, int]],
    ads_edges: set[tuple[int, int]],
    host_edges: set[tuple[int, int]],
    *,
    center_elements: tuple[str, ...],
    framework_elements: tuple[str, ...],
    id_prefix: str,
    family: str,
    per_center_naming: bool,
    allowed_state_elements: tuple[str, ...],
    max_state_members: int | None,
    rule_profile: str,
    output_index_base: int,
    start_index: int = 1,
) -> list[GeneratedActiveSite]:
    syms = atoms.get_chemical_symbols()
    cov_adj = _adjacency(cov_edges)
    ads_adj = _adjacency(ads_edges)
    host_adj = _adjacency(host_edges)
    structural_adj = _adjacency(cov_edges | ads_edges | host_edges)

    sites: list[GeneratedActiveSite] = []
    next_idx = int(start_index)
    next_idx_by_element: dict[str, int] = defaultdict(lambda: int(start_index))

    center_set = {x.strip() for x in center_elements if x.strip()}
    framework_set = {x.strip() for x in framework_elements if x.strip()}

    for center in range(len(atoms)):
        if not bool(host_mask[center]):
            continue
        if syms[center] not in center_set:
            continue

        oxygen_candidates = sorted(
            nbr
            for nbr in structural_adj.get(center, set())
            if syms[nbr] == "O"
        )
        for oxygen in oxygen_candidates:
            h_neighbors = sorted(n for n in structural_adj.get(oxygen, set()) if syms[n] == "H")
            if not h_neighbors:
                continue

            host_neighbors = sorted(n for n in host_adj.get(oxygen, set()) if bool(host_mask[n]))
            framework_candidates = [n for n in host_neighbors if n != center and syms[n] in framework_set]
            if not framework_candidates:
                continue

            framework_atom = _pick_best_framework_neighbor(atoms, oxygen, framework_candidates)

            for hydrogen in h_neighbors:
                center_element = syms[center]
                if per_center_naming:
                    local_prefix = f"{id_prefix}_{center_element}"
                    local_family = f"{family}_{center_element}"
                    seq = next_idx_by_element[center_element]
                    site_id = f"{local_prefix}_{seq:02d}"
                    next_idx_by_element[center_element] = seq + 1
                else:
                    local_prefix = id_prefix
                    local_family = family
                    site_id = f"{local_prefix}_{next_idx:02d}"
                    next_idx += 1
                sites.append(
                    GeneratedActiveSite(
                        site_id=site_id,
                        family=local_family,
                        core_members=(center, oxygen, framework_atom),
                        initial_state_members=(hydrogen,),
                        allowed_state_elements=tuple(allowed_state_elements),
                        max_state_members=max_state_members,
                        rule_profile=rule_profile,
                        metadata={
                            "generator": "bas_oh_bridge",
                            "center_element": syms[center],
                            "center_atom": int(center) + (1 if output_index_base == 1 else 0),
                            "oxygen_atom": int(oxygen) + (1 if output_index_base == 1 else 0),
                            "framework_atom": int(framework_atom) + (1 if output_index_base == 1 else 0),
                            "framework_element": syms[framework_atom],
                            "initial_hydrogen": int(hydrogen) + (1 if output_index_base == 1 else 0),
                        },
                    )
                )

    return sites


def _dump_yaml(path: str, data: dict[str, Any]) -> None:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError("Writing active-site YAML requires PyYAML. Please `pip install pyyaml`.") from e

    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def _default_output_path(traj_path: str) -> str:
    p = Path(traj_path)
    return str(p.with_name(f"{p.stem}_active_site.yaml"))


def build_active_sites(
    *,
    config_path: str,
    structure_path: str | None,
    traj_override: str | None,
    frame_index: int,
    motif: str,
    center_elements: list[str],
    framework_elements: list[str],
    id_prefix: str,
    family: str,
    per_center_naming: bool,
    allowed_state_elements: list[str],
    max_state_members: int | None,
    rule_profile: str,
    output_path: str | None,
    output_index_base: int,
    start_index: int,
) -> str:
    cfg, inferred_base = load_raw_config(config_path)
    merged = apply_overrides(cfg, RunOverrides(traj=traj_override))
    criteria = _build_criteria(merged)
    host_def = _build_slab_def(merged)

    if structure_path:
        structure_source = _resolve_path(structure_path, inferred_base)
    else:
        run_cfg = _get_section(merged, "run")
        traj_raw = str(run_cfg.get("traj", "")).strip()
        if not traj_raw:
            raise ValueError(
                "Missing STRUCTURE_PATH and run.traj in config. "
                "Please provide a first-frame structure file or configure run.traj."
            )
        structure_source = _resolve_path(traj_raw, inferred_base)

    atoms = _load_structure(structure_source, frame_index)
    host_mask = host_def.mask(atoms)
    cov_edges, ads_edges, host_edges = _build_edges(atoms, host_mask, criteria)

    motif_key = str(motif).strip().lower()
    if motif_key not in {"bas_oh_bridge", "bas"}:
        raise ValueError(f"Unsupported MOTIF: {motif}. Currently supported: bas_oh_bridge")

    sites = generate_bas_oh_bridge_sites(
        atoms,
        host_mask,
        cov_edges,
        ads_edges,
        host_edges,
        center_elements=tuple(center_elements),
        framework_elements=tuple(framework_elements),
        id_prefix=id_prefix,
        family=family,
        per_center_naming=per_center_naming,
        allowed_state_elements=tuple(allowed_state_elements),
        max_state_members=max_state_members,
        rule_profile=rule_profile,
        output_index_base=output_index_base,
        start_index=start_index,
    )

    out_path = output_path or _default_output_path(structure_source)
    payload = {"sites": [s.to_mapping(output_index_base=output_index_base) for s in sites]}
    _dump_yaml(out_path, payload)
    print(
        f"[info] structure={structure_source} host_atoms={int(sum(bool(x) for x in host_mask))} "
        f"sites_found={len(sites)}"
    )
    if not sites:
        print(
            "[warn] No active sites were generated. Check whether:\n"
            "  1. host_definition marks center/framework atoms as host,\n"
            "  2. the motif really looks like center-O(H)-framework in the provided structure,\n"
            "  3. CENTER_ELEMENTS / FRAMEWORK_ELEMENTS match the actual elements."
        )
    return out_path


if __name__ == "__main__":
    # --------------------------------------------
    # Edit these parameters directly in your IDE.
    # --------------------------------------------
    # CONFIG_PATH is still used to read host_definition and criteria.
    CONFIG_PATH = "GaH2_gfn_xtb_analyzer.yaml"

    # Preferred: point this to a first-frame structure file so we do not have to
    # load the full trajectory just to generate active-site seeds.
    # Examples: .xyz, .extxyz, POSCAR, CONTCAR, .traj
    STRUCTURE_PATH = 'example_traj/GaH2_gfn_xtb/GaH2_gfn_xtb_0.cif'

    # Optional fallback: override run.traj in CONFIG_PATH when STRUCTURE_PATH is None.
    TRAJ_OVERRIDE = None

    # FRAME_INDEX only matters for multi-frame files.
    FRAME_INDEX = 0

    # Currently supported:
    # - "bas_oh_bridge": center - O(H) - framework
    MOTIF = "bas_oh_bridge"

    # Multiple center elements are allowed, for example:
    # ["Ga", "Al"]
    CENTER_ELEMENTS = ["Ga", 'Al']

    # Framework elements bonded to the acidic oxygen.
    FRAMEWORK_ELEMENTS = ["Si"]

    ID_PREFIX = "B"
    FAMILY = "BAS"
    # If True, generated ids/families are automatically split by center element.
    # Example:
    # CENTER_ELEMENTS = ["Ga", "Al"]
    # ID_PREFIX = "B"
    # FAMILY = "BAS"
    # -> B_Ga_01 / BAS_Ga, B_Al_01 / BAS_Al
    PER_CENTER_NAMING = True
    ALLOWED_STATE_ELEMENTS = ["H", "O", "C"]
    MAX_STATE_MEMBERS = 20
    RULE_PROFILE = "bas_default"

    # Output YAML index style:
    # 0 = 0-based
    # 1 = 1-based
    OUTPUT_INDEX_BASE = 1
    START_INDEX = 1
    OUTPUT_PATH = "BAS_active_site.yaml"

    out = build_active_sites(
        config_path=CONFIG_PATH,
        structure_path=STRUCTURE_PATH,
        traj_override=TRAJ_OVERRIDE,
        frame_index=FRAME_INDEX,
        motif=MOTIF,
        center_elements=CENTER_ELEMENTS,
        framework_elements=FRAMEWORK_ELEMENTS,
        id_prefix=ID_PREFIX,
        family=FAMILY,
        per_center_naming=PER_CENTER_NAMING,
        allowed_state_elements=ALLOWED_STATE_ELEMENTS,
        max_state_members=MAX_STATE_MEMBERS,
        rule_profile=RULE_PROFILE,
        output_path=OUTPUT_PATH,
        output_index_base=OUTPUT_INDEX_BASE,
        start_index=START_INDEX,
    )
    print(f"[done] active-site YAML written to: {out}")
