"""Story-mode source-to-target path extraction."""

from __future__ import annotations

import math
from typing import Any

import networkx as nx

from .base import (
    collapse_reversible,
    filter_edges_by_min_weight,
    node_name_to_id,
    prune_isolates,
    summarize_graph,
)


def _resolve_nodes(graph: nx.Graph, raw_names: list[str], *, role: str) -> list[str]:
    out: list[str] = []
    for name in raw_names:
        node_id = node_name_to_id(graph, name)
        if node_id is None:
            print(f"[story] WARNING: {role} not found: {name}")
            continue
        out.append(node_id)
    return out


def _species_formula(node_data: dict[str, Any]) -> str:
    return str(node_data.get("species_formula", "") or "").strip()


def _species_label(node_data: dict[str, Any]) -> str:
    return str(node_data.get("orig_id", "") or node_data.get("species_label", "") or "").strip()


def _filter_story_intermediates(
    graph: nx.Graph,
    *,
    protected_nodes: set[str],
    excluded_species_formulas: set[str],
    excluded_species_labels: set[str],
) -> nx.Graph:
    if not excluded_species_formulas and not excluded_species_labels:
        return graph.copy()

    out = graph.copy()
    to_remove: list[str] = []
    for node_id, data in out.nodes(data=True):
        node_id_str = str(node_id)
        if node_id_str in protected_nodes:
            continue
        if str(data.get("node_type", "")) != "species":
            continue
        formula = _species_formula(data)
        label = _species_label(data)
        if formula in excluded_species_formulas or label in excluded_species_labels:
            to_remove.append(node_id_str)
    if to_remove:
        out.remove_nodes_from(to_remove)
    return out


def _build_story_work_graph(
    graph: nx.Graph,
    *,
    direction: str,
    path_mode: str,
    hub_penalty_strength: float,
    source_nodes: set[str],
    target_nodes: set[str],
) -> nx.Graph:
    direction_norm = str(direction or "directed").strip().lower()
    path_mode_norm = str(path_mode or "shortest").strip().lower()
    work_graph = graph if (graph.is_directed() and direction_norm == "directed") else graph.to_undirected()

    if path_mode_norm == "shortest":
        return work_graph
    if path_mode_norm != "chemical_shortest":
        raise ValueError("story mode currently supports path_mode='shortest' or 'chemical_shortest'.")

    weighted = work_graph.copy()
    source_target = set(source_nodes) | set(target_nodes)

    for u, v, data in weighted.edges(data=True):
        cost = 1.0
        target_node = v
        target_data = weighted.nodes[target_node]
        if str(target_data.get("node_type", "")) == "species" and str(target_node) not in source_target:
            degree = weighted.degree(target_node)
            # Penalize hub-like species so that ubiquitous shuttles such as H are less likely
            # to dominate the extracted story unless no better path exists.
            cost += max(0.0, math.log2(max(1, degree)) - 1.0) * float(hub_penalty_strength)
        data["story_cost"] = cost

    return weighted


def run_story_mode(graph: nx.Graph, cfg: dict[str, Any]) -> nx.Graph:
    story_cfg = cfg.get("story", {})
    out = graph.copy()

    min_edge_weight = int(story_cfg.get("min_edge_weight", 1) or 1)
    out = filter_edges_by_min_weight(out, min_edge_weight)
    summarize_graph(out, "story:min_edge_weight")

    if bool(story_cfg.get("collapse_reversible", False)):
        out = collapse_reversible(out)
        summarize_graph(out, "story:collapse_reversible")

    sources = _resolve_nodes(out, list(story_cfg.get("sources", []) or []), role="source")
    targets = _resolve_nodes(out, list(story_cfg.get("targets", []) or []), role="target")
    if not sources or not targets:
        raise ValueError("story mode requires at least one resolvable source and one resolvable target.")

    direction = str(story_cfg.get("direction", "directed") or "directed").strip().lower()
    path_mode = str(story_cfg.get("path_mode", "shortest") or "shortest").strip().lower()
    max_paths = int(story_cfg.get("max_paths", 10) or 10)
    if max_paths <= 0:
        raise ValueError("story.max_paths must be >= 1.")

    excluded_species_formulas = {
        str(x).strip() for x in list(story_cfg.get("excluded_species_formulas", []) or []) if str(x).strip()
    }
    excluded_species_labels = {
        str(x).strip() for x in list(story_cfg.get("excluded_species_labels", []) or []) if str(x).strip()
    }
    out = _filter_story_intermediates(
        out,
        protected_nodes=set(sources) | set(targets),
        excluded_species_formulas=excluded_species_formulas,
        excluded_species_labels=excluded_species_labels,
    )
    summarize_graph(out, "story:exclude_intermediates")

    hub_penalty_strength = float(story_cfg.get("hub_penalty_strength", 2.0) or 0.0)
    work_graph = _build_story_work_graph(
        out,
        direction=direction,
        path_mode=path_mode,
        hub_penalty_strength=hub_penalty_strength,
        source_nodes=set(sources),
        target_nodes=set(targets),
    )

    keep_nodes: set[str] = set()
    keep_edges: set[tuple[str, str]] = set()
    kept_paths = 0

    for source in sources:
        for target in targets:
            if kept_paths >= max_paths:
                break
            try:
                kwargs: dict[str, Any] = {}
                if path_mode == "chemical_shortest":
                    kwargs["weight"] = "story_cost"
                for path in nx.all_shortest_paths(work_graph, source=source, target=target, **kwargs):
                    keep_nodes.update(str(node_id) for node_id in path)
                    keep_edges.update((str(u), str(v)) for u, v in zip(path[:-1], path[1:]))
                    kept_paths += 1
                    if kept_paths >= max_paths:
                        break
            except nx.NetworkXNoPath:
                continue
        if kept_paths >= max_paths:
            break

    if not keep_nodes:
        raise ValueError("story mode found no path between the requested sources and targets.")

    if out.is_directed():
        directed_edges: set[tuple[str, str]] = set()
        for u, v in keep_edges:
            if out.has_edge(u, v):
                directed_edges.add((u, v))
            if direction != "directed" and out.has_edge(v, u):
                directed_edges.add((v, u))
        subgraph = out.edge_subgraph(directed_edges).copy() if directed_edges else out.subgraph(keep_nodes).copy()
        if subgraph.number_of_nodes() == 0:
            subgraph = out.subgraph(keep_nodes).copy()
    else:
        subgraph = out.subgraph(keep_nodes).copy()

    summarize_graph(subgraph, "story:selected")

    if bool(story_cfg.get("prune_isolates", True)):
        subgraph = prune_isolates(subgraph)
        summarize_graph(subgraph, "story:prune_isolates")

    return subgraph
