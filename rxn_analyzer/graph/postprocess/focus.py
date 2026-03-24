"""Focus-mode graph simplification."""

from __future__ import annotations

from typing import Any

import networkx as nx

from .base import (
    collapse_reversible,
    compute_focus_scores,
    filter_edges_by_min_weight,
    node_name_to_id,
    prune_isolates,
    summarize_graph,
    write_focus_scores_csv,
)


def _rank_nodes(graph: nx.Graph, scores: dict[str, int], node_type: str) -> list[str]:
    ranked: list[tuple[int, str]] = []
    for node_id, data in graph.nodes(data=True):
        if str(data.get("node_type", "")) != node_type:
            continue
        ranked.append((scores.get(str(node_id), 0), str(node_id)))
    ranked.sort(key=lambda item: (-item[0], item[1]))
    return [node_id for _score, node_id in ranked]


def _species_family(node_data: dict[str, Any], mode: str) -> str:
    mode_norm = str(mode or "ads_state").strip().lower()
    if mode_norm == "ads_state":
        return "ads" if bool(node_data.get("species_is_ads", node_data.get("is_ads", False))) else "non_ads"
    raise ValueError(f"Unknown species_family_mode: {mode}. Supported modes: ads_state.")


def _select_species_with_family_quota(
    graph: nx.Graph,
    ranked_species: list[str],
    *,
    top_species: int,
    top_species_per_family: int,
    species_family_mode: str,
) -> set[str]:
    if top_species <= 0:
        return set()
    if top_species_per_family <= 0:
        return set(ranked_species[:top_species])

    selected: list[str] = []
    selected_set: set[str] = set()
    family_counts: dict[str, int] = {}

    for node_id in ranked_species:
        family = _species_family(graph.nodes[node_id], species_family_mode)
        if family_counts.get(family, 0) >= top_species_per_family:
            continue
        selected.append(node_id)
        selected_set.add(node_id)
        family_counts[family] = family_counts.get(family, 0) + 1

    if len(selected) < top_species:
        for node_id in ranked_species:
            if node_id in selected_set:
                continue
            selected.append(node_id)
            selected_set.add(node_id)
            if len(selected) >= top_species:
                break
    else:
        selected = selected[:top_species]
        selected_set = set(selected)

    return selected_set


def _resolve_seed_nodes(graph: nx.Graph, raw_seeds: list[str]) -> set[str]:
    out: set[str] = set()
    for seed in raw_seeds:
        node_id = node_name_to_id(graph, seed)
        if node_id is None:
            print(f"[focus] WARNING: seed not found: {seed}")
            continue
        out.add(node_id)
    return out


def _expand_neighbors(graph: nx.Graph, seeds: set[str], depth: int) -> set[str]:
    if depth <= 0 or not seeds:
        return set(seeds)
    visited = set(seeds)
    frontier = set(seeds)
    for _ in range(depth):
        nxt: set[str] = set()
        for node_id in frontier:
            nxt.update(str(nei) for nei in graph.neighbors(node_id))
            if graph.is_directed():
                nxt.update(str(nei) for nei in graph.predecessors(node_id))
        nxt -= visited
        visited |= nxt
        frontier = nxt
        if not frontier:
            break
    return visited


def run_focus_mode(
    graph: nx.Graph,
    cfg: dict[str, Any],
    *,
    score_csv_path: str | None = None,
) -> nx.Graph:
    focus_cfg = cfg.get("focus", {})
    out = graph.copy()

    min_edge_weight = int(focus_cfg.get("min_edge_weight", 1) or 1)
    out = filter_edges_by_min_weight(out, min_edge_weight)
    summarize_graph(out, "focus:min_edge_weight")

    if bool(focus_cfg.get("collapse_reversible", False)):
        out = collapse_reversible(out)
        summarize_graph(out, "focus:collapse_reversible")

    score_mode = str(focus_cfg.get("score_mode", "weighted_degree") or "weighted_degree")
    scores = compute_focus_scores(out, score_mode)
    ranked_species = _rank_nodes(out, scores, "species")
    ranked_reactions = _rank_nodes(out, scores, "reaction")

    keep_species = int(focus_cfg.get("top_species", 30) or 30)
    keep_reactions = int(focus_cfg.get("top_reactions", 20) or 20)
    top_species_per_family = int(focus_cfg.get("top_species_per_family", 0) or 0)
    species_family_mode = str(focus_cfg.get("species_family_mode", "ads_state") or "ads_state")
    include_neighbors = bool(focus_cfg.get("include_neighbors", True))
    prune = bool(focus_cfg.get("prune_isolates", True))
    seed_nodes = _resolve_seed_nodes(out, list(focus_cfg.get("keep_nodes", []) or []))
    protect_seed_neighbors = int(focus_cfg.get("protect_seed_neighbors", 0) or 0)

    selected_species = _select_species_with_family_quota(
        out,
        ranked_species,
        top_species=keep_species,
        top_species_per_family=top_species_per_family,
        species_family_mode=species_family_mode,
    )
    selected: set[str] = selected_species | set(ranked_reactions[:keep_reactions]) | seed_nodes

    if protect_seed_neighbors > 0:
        selected |= _expand_neighbors(out, seed_nodes, protect_seed_neighbors)

    if include_neighbors:
        selected = _expand_neighbors(out, selected, 1)

    focused = out.subgraph(selected).copy()
    summarize_graph(focused, "focus:selected")

    if prune:
        focused = prune_isolates(focused)
        summarize_graph(focused, "focus:prune_isolates")

    if score_csv_path:
        write_focus_scores_csv(
            score_csv_path,
            out,
            scores,
            score_mode=score_mode,
            selected=set(focused.nodes()),
        )

    return focused
