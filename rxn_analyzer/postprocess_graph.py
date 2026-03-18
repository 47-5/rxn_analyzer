#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from typing import Dict, Set, Tuple, List

import networkx as nx


# -----------------------
# Helpers
# -----------------------
def safe_int(v, default=0) -> int:
    try:
        return int(v)
    except Exception:
        return default


def edge_types(d: dict) -> Set[str]:
    t = d.get("type", "")
    if not t:
        return set()
    return {x.strip() for x in str(t).split("|") if x.strip()}


def node_name_to_id(G: nx.Graph, name: str) -> str | None:
    # Try direct node id
    if G.has_node(name):
        return str(name)

    # Try by orig_id
    for n, d in G.nodes(data=True):
        if str(d.get("orig_id")) == name:
            return str(n)

    # Try by label
    for n, d in G.nodes(data=True):
        if str(d.get("label")) == name:
            return str(n)

    return None


def labelize(long_id: str) -> str:
    return long_id.replace(" + ", " +\n").replace("|ads", "\n|ads")


def ensure_node(G: nx.Graph, long_id: str) -> str:
    # If node exists by orig_id, return it
    for n, d in G.nodes(data=True):
        if str(d.get("orig_id")) == long_id:
            return str(n)

    # Otherwise create a new node with minimal attrs
    new_id = str(len(G.nodes))
    G.add_node(new_id)
    G.nodes[new_id]["orig_id"] = long_id
    G.nodes[new_id]["label"] = labelize(long_id)
    return new_id


def summarize(G: nx.Graph, tag: str):
    print(f"[{tag}] nodes={G.number_of_nodes()} edges={G.number_of_edges()}")


# -----------------------
# Event aggregation (for time window / node freq)
# -----------------------
def load_events_transform_csv(path: str) -> List[dict]:
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)
    return rows


def aggregate_events(
    events: List[dict],
    frame_min: int | None,
    frame_max: int | None,
    type_filter: Set[str] | None = None,
) -> tuple[Dict[Tuple[str, str], int], Dict[Tuple[str, str], Set[str]], Dict[str, int]]:
    """
    Aggregate TransformEvents into directed edge counts.
    Returns:
      edge_counts[(src_long, dst_long)] = count
      edge_types[(src_long, dst_long)] = {type,...}
      node_freq[long_id] = frequency in reactants+products
    """
    edge_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    edge_types: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    node_freq: Dict[str, int] = defaultdict(int)

    for row in events:
        frame = safe_int(row.get("frame", -1))
        if frame_min is not None and frame < frame_min:
            continue
        if frame_max is not None and frame > frame_max:
            continue

        ttype = str(row.get("type", "") or "")
        if type_filter is not None and ttype not in type_filter:
            continue

        reactants = json.loads(row.get("reactants", "[]"))
        products = json.loads(row.get("products", "[]"))

        reactants_sorted = sorted(reactants)
        products_sorted = sorted(products)

        long_src = " + ".join(reactants_sorted)
        long_dst = " + ".join(products_sorted)

        edge_counts[(long_src, long_dst)] += 1
        if ttype:
            edge_types[(long_src, long_dst)].add(ttype)

        for r in reactants_sorted:
            node_freq[r] += 1
        for p in products_sorted:
            node_freq[p] += 1

    return edge_counts, edge_types, node_freq


# -----------------------
# Pipeline Ops
# -----------------------
def op_time_window(G: nx.DiGraph, cfg: dict, events_csv: str | None, node_freq_holder: dict) -> nx.DiGraph:
    if not events_csv:
        print("[time_window] WARNING: events CSV not provided, skipping time window.")
        return G

    frame_min = cfg.get("start", None)
    frame_max = cfg.get("end", None)
    type_filter = set(cfg.get("keep_types", [])) or None

    events = load_events_transform_csv(events_csv)
    edge_counts, edge_types, node_freq = aggregate_events(
        events, frame_min=frame_min, frame_max=frame_max, type_filter=type_filter
    )

    node_freq_holder.clear()
    node_freq_holder.update(node_freq)

    # rebuild edge set based on window counts
    H = G.copy()
    H.remove_edges_from(list(H.edges()))

    # add edges from aggregated events
    for (src_long, dst_long), w in edge_counts.items():
        if w <= 0:
            continue
        src = ensure_node(H, src_long)
        dst = ensure_node(H, dst_long)
        types = edge_types.get((src_long, dst_long), set())
        ttype = "|".join(sorted(types)) if types else ""
        H.add_edge(src, dst, weight=int(w), type=ttype, label=f"{ttype} ({w})")

    summarize(H, "time_window")
    return H


def op_filter_type(G: nx.DiGraph, cfg: dict) -> nx.DiGraph:
    keep = set(cfg.get("keep", []))
    if not keep:
        return G

    H = G.copy()
    for u, v, d in list(H.edges(data=True)):
        types = edge_types(d)
        if not types.intersection(keep):
            H.remove_edge(u, v)

    summarize(H, "filter_type")
    return H


def op_min_edge_weight(G: nx.DiGraph, cfg: dict) -> nx.DiGraph:
    m = safe_int(cfg.get("min", 1), 1)
    H = G.copy()
    for u, v, d in list(H.edges(data=True)):
        w = safe_int(d.get("weight", 1), 1)
        if w < m:
            H.remove_edge(u, v)
    summarize(H, "min_edge_weight")
    return H


def op_min_node_freq(G: nx.DiGraph, cfg: dict, node_freq_holder: dict) -> nx.DiGraph:
    m = safe_int(cfg.get("min", 1), 1)
    source = cfg.get("source", "events")  # events | degree
    use_degree_proxy = bool(cfg.get("use_degree_proxy", True))

    freq: Dict[str, int] = {}

    if source == "events" and node_freq_holder:
        # node_freq_holder keys are long_id; map to node ids
        for n, d in G.nodes(data=True):
            long_id = str(d.get("orig_id", ""))
            if long_id in node_freq_holder:
                freq[n] = node_freq_holder[long_id]
            else:
                freq[n] = 0
    else:
        if not use_degree_proxy:
            print("[min_node_freq] WARNING: no event freq and proxy disabled; skipping.")
            return G
        # weighted degree proxy
        for n in G.nodes():
            deg = 0
            for _u, _v, d in G.in_edges(n, data=True):
                deg += safe_int(d.get("weight", 1), 1)
            for _u, _v, d in G.out_edges(n, data=True):
                deg += safe_int(d.get("weight", 1), 1)
            freq[n] = deg

    H = G.copy()
    for n in list(H.nodes()):
        if freq.get(n, 0) < m:
            H.remove_node(n)

    summarize(H, "min_node_freq")
    return H


def op_k_core(G: nx.DiGraph, cfg: dict) -> nx.DiGraph:
    k = safe_int(cfg.get("k", 2), 2)
    mode = cfg.get("mode", "undirected")  # undirected | directed
    if mode == "undirected":
        U = G.to_undirected()
        core = nx.k_core(U, k=k)
        H = G.subgraph(core.nodes()).copy()
    else:
        # directed k-core not standard; use undirected fallback
        U = G.to_undirected()
        core = nx.k_core(U, k=k)
        H = G.subgraph(core.nodes()).copy()

    summarize(H, "k_core")
    return H


def op_neighborhood(G: nx.DiGraph, cfg: dict) -> nx.DiGraph:
    seeds = cfg.get("seeds", [])
    depth = safe_int(cfg.get("depth", 1), 1)
    mode = cfg.get("mode", "undirected")  # undirected | out | in

    seed_ids = []
    for s in seeds:
        nid = node_name_to_id(G, s)
        if nid is not None:
            seed_ids.append(nid)
        else:
            print(f"[neighborhood] WARNING: seed not found: {s}")

    if not seed_ids:
        return G

    if mode == "undirected":
        H = G.to_undirected()
        frontier = set(seed_ids)
        visited = set(seed_ids)
        for _ in range(depth):
            nxt = set()
            for u in frontier:
                nxt.update(H.neighbors(u))
            nxt -= visited
            visited |= nxt
            frontier = nxt
        sub = G.subgraph(visited).copy()
    elif mode == "out":
        frontier = set(seed_ids)
        visited = set(seed_ids)
        for _ in range(depth):
            nxt = set()
            for u in frontier:
                nxt.update(v for _, v in G.out_edges(u))
            nxt -= visited
            visited |= nxt
            frontier = nxt
        sub = G.subgraph(visited).copy()
    else:  # in
        frontier = set(seed_ids)
        visited = set(seed_ids)
        for _ in range(depth):
            nxt = set()
            for u in frontier:
                nxt.update(v for v, _ in G.in_edges(u))
            nxt -= visited
            visited |= nxt
            frontier = nxt
        sub = G.subgraph(visited).copy()

    summarize(sub, "neighborhood")
    return sub


def op_shortest_paths(G: nx.DiGraph, cfg: dict) -> nx.DiGraph:
    sources = cfg.get("sources", [])
    targets = cfg.get("targets", [])
    mode = cfg.get("mode", "directed")  # directed | undirected

    src_ids = []
    tgt_ids = []
    for s in sources:
        nid = node_name_to_id(G, s)
        if nid is not None:
            src_ids.append(nid)
        else:
            print(f"[shortest_paths] WARNING: source not found: {s}")
    for t in targets:
        nid = node_name_to_id(G, t)
        if nid is not None:
            tgt_ids.append(nid)
        else:
            print(f"[shortest_paths] WARNING: target not found: {t}")

    if not src_ids or not tgt_ids:
        return G

    H = G if mode == "directed" else G.to_undirected()

    keep_nodes: Set[str] = set()
    keep_edges: Set[Tuple[str, str]] = set()

    for s in src_ids:
        for t in tgt_ids:
            try:
                paths = list(nx.all_shortest_paths(H, source=s, target=t))
                for p in paths:
                    keep_nodes.update(p)
                    keep_edges.update(zip(p[:-1], p[1:]))
            except nx.NetworkXNoPath:
                continue

    if mode == "undirected":
        sub = H.edge_subgraph(keep_edges).copy()
        sub = G.subgraph(sub.nodes()).copy()
    else:
        sub = G.edge_subgraph(keep_edges).copy()

    summarize(sub, "shortest_paths")
    return sub


def op_collapse_reversible(G: nx.DiGraph, cfg: dict) -> nx.Graph:
    U = nx.Graph()
    # copy nodes
    for n, d in G.nodes(data=True):
        U.add_node(n, **d)

    for u, v, d in G.edges(data=True):
        w = safe_int(d.get("weight", 1), 1)
        types = edge_types(d)
        if U.has_edge(u, v):
            U[u][v]["weight"] = safe_int(U[u][v].get("weight", 0), 0) + w
            old_types = set(str(U[u][v].get("type", "")).split("|")) if U[u][v].get("type") else set()
            new_types = old_types | types
            U[u][v]["type"] = "|".join(sorted(t for t in new_types if t))
        else:
            U.add_edge(u, v, weight=w, type="|".join(sorted(types)))

    # update labels
    for u, v, d in U.edges(data=True):
        d["label"] = f"{d.get('type','')} ({d.get('weight',1)})"

    summarize(U, "collapse_reversible")
    return U


def prune_isolates(G: nx.Graph) -> nx.Graph:
    H = G.copy()
    isolates = list(nx.isolates(H))
    H.remove_nodes_from(isolates)
    summarize(H, "prune_isolates")
    return H


# -----------------------
# Public API (package use)
# -----------------------
def load_config(cfg: dict | str) -> dict:
    """
    Accepts:
      - dict (already parsed)
      - str  (path to JSON)
    """
    if isinstance(cfg, dict):
        return cfg
    with open(cfg, "r", encoding="utf-8") as f:
        return json.load(f)


def run_pipeline(
    G: nx.Graph,
    cfg: dict | str,
    *,
    events_csv: str | None = None,
    prune_isolates_flag: bool | None = None,
    verbose: bool = True,
) -> nx.Graph:
    """
    Run the configurable postprocess pipeline on a graph.

    Parameters
    ----------
    G : nx.Graph
        Input graph (usually DiGraph from GraphML).
    cfg : dict | str
        Config dict or JSON path.
    events_csv : str | None
        Path to transform events CSV (for time window / node freq).
    prune_isolates_flag : bool | None
        Override config's prune_isolates if not None.
    verbose : bool
        Print summary info.

    Returns
    -------
    nx.Graph
        Processed graph (may be DiGraph or Graph depending on pipeline).
    """
    cfg = load_config(cfg)
    if verbose:
        summarize(G, "input")

    node_freq_holder: Dict[str, int] = {}

    for step in cfg.get("pipeline", []):
        op = step.get("op", "")
        if op == "time_window":
            G = op_time_window(G, step, events_csv, node_freq_holder)
        elif op == "filter_type":
            G = op_filter_type(G, step)
        elif op == "min_edge_weight":
            G = op_min_edge_weight(G, step)
        elif op == "min_node_freq":
            G = op_min_node_freq(G, step, node_freq_holder)
        elif op == "k_core":
            G = op_k_core(G, step)
        elif op == "neighborhood":
            G = op_neighborhood(G, step)
        elif op == "shortest_paths":
            G = op_shortest_paths(G, step)
        elif op == "collapse_reversible":
            G = op_collapse_reversible(G, step)
        else:
            if verbose:
                print(f"[pipeline] WARNING: unknown op: {op}")

    prune_flag = cfg.get("prune_isolates", True) if prune_isolates_flag is None else prune_isolates_flag
    if bool(prune_flag):
        G = prune_isolates(G)

    return G


# -----------------------
# CLI
# -----------------------
def main():
    ap = argparse.ArgumentParser("postprocess reaction network graph")
    ap.add_argument("--in-graph", required=True, help="Input GraphML")
    ap.add_argument("--out-graph", required=True, help="Output GraphML")
    ap.add_argument("--config", required=True, help="Config JSON")
    ap.add_argument("--events", default=None, help="Transform events CSV (for time window / node freq)")
    args = ap.parse_args()

    G = nx.read_graphml(args.in_graph)
    H = run_pipeline(G, args.config, events_csv=args.events, verbose=True)

    nx.write_graphml(H, args.out_graph)
    print(f"[done] wrote {args.out_graph}")


if __name__ == "__main__":
    main()