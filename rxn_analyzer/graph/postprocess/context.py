"""Context-mode local subgraph extraction."""

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


def _resolve_seed_nodes(graph: nx.Graph, raw_seeds: list[str]) -> set[str]:
    out: set[str] = set()
    for seed in raw_seeds:
        node_id = node_name_to_id(graph, seed)
        if node_id is None:
            print(f"[context] WARNING: seed not found: {seed}")
            continue
        out.add(node_id)
    return out


def _expand_context(graph: nx.Graph, seeds: set[str], depth: int, direction: str) -> set[str]:
    if depth <= 0 or not seeds:
        return set(seeds)

    direction_norm = str(direction or "both").strip().lower()
    visited = set(seeds)
    frontier = set(seeds)

    for _ in range(depth):
        nxt: set[str] = set()
        for node_id in frontier:
            if direction_norm in {"both", "out"}:
                nxt.update(str(nei) for nei in graph.neighbors(node_id))
            if graph.is_directed() and direction_norm in {"both", "in"}:
                nxt.update(str(nei) for nei in graph.predecessors(node_id))
            if not graph.is_directed() and direction_norm == "in":
                nxt.update(str(nei) for nei in graph.neighbors(node_id))
        nxt -= visited
        visited |= nxt
        frontier = nxt
        if not frontier:
            break

    return visited


def run_context_mode(graph: nx.Graph, cfg: dict[str, Any]) -> nx.Graph:
    context_cfg = cfg.get("context", {})
    out = graph.copy()

    min_edge_weight = int(context_cfg.get("min_edge_weight", 1) or 1)
    out = filter_edges_by_min_weight(out, min_edge_weight)
    summarize_graph(out, "context:min_edge_weight")

    if bool(context_cfg.get("collapse_reversible", False)):
        out = collapse_reversible(out)
        summarize_graph(out, "context:collapse_reversible")

    seeds = _resolve_seed_nodes(out, list(context_cfg.get("seeds", []) or []))
    if not seeds:
        raise ValueError("context mode requires at least one resolvable seed in context.seeds.")

    depth = int(context_cfg.get("depth", 1) or 1)
    direction = str(context_cfg.get("direction", "both") or "both")
    selected = _expand_context(out, seeds, depth, direction)

    subgraph = out.subgraph(selected).copy()
    summarize_graph(subgraph, "context:selected")

    if bool(context_cfg.get("prune_isolates", True)):
        subgraph = prune_isolates(subgraph)
        summarize_graph(subgraph, "context:prune_isolates")

    return subgraph
