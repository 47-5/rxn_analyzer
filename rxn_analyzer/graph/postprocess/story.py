"""Story-mode source-to-target path extraction."""

from __future__ import annotations

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
    if path_mode != "shortest":
        raise ValueError("story mode currently supports only path_mode='shortest'.")

    max_paths = int(story_cfg.get("max_paths", 10) or 10)
    if max_paths <= 0:
        raise ValueError("story.max_paths must be >= 1.")

    work_graph = out if (out.is_directed() and direction == "directed") else out.to_undirected()

    keep_nodes: set[str] = set()
    keep_edges: set[tuple[str, str]] = set()
    kept_paths = 0

    for source in sources:
        for target in targets:
            if kept_paths >= max_paths:
                break
            try:
                for path in nx.all_shortest_paths(work_graph, source=source, target=target):
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
