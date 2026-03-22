"""Main reaction-graph construction helpers.

This module is intentionally limited to the ordinary reaction-network graph:
- create species and reaction nodes
- add transform events to the graph
- summarize reversible reaction pairs

Low-level ID allocation lives in `graph.ids`.
Attribute parsing/formatting lives in `graph.attrs`.
"""

from __future__ import annotations

import json
from collections import Counter

import networkx as nx

from .attrs import (
    build_label_fields,
    node_attrs_from_species_list,
    parse_frames_attr,
    parse_species_label,
)
from .ids import (
    get_or_create_numeric_node_id,
    get_or_create_reaction_node_id,
    reaction_signature,
)


def ensure_graph() -> nx.DiGraph:
    return nx.DiGraph()


def ensure_bipartite_graph() -> nx.DiGraph:
    return nx.DiGraph()


def ensure_species_node(graph: nx.DiGraph, long_id: str, species: list[str]) -> str:
    node_id = get_or_create_numeric_node_id(long_id)
    if not graph.has_node(node_id):
        graph.add_node(node_id)

    data = graph.nodes[node_id]
    parsed = [parse_species_label(label) for label in species]
    short_display = " + ".join(str(part["short"]) for part in parsed) if parsed else long_id

    data.setdefault("orig_id", long_id)
    data.setdefault("label", long_id.replace(" + ", " +\n").replace("|ads", "\n|ads"))
    data.setdefault("display_label", short_display)
    data.setdefault("node_type", "species")
    data.setdefault("is_reaction", False)
    data.setdefault("bipartite", "species")
    data.setdefault("species_labels", json.dumps(species, ensure_ascii=False))

    for key, value in node_attrs_from_species_list(species).items():
        data.setdefault(key, value)
    for key, value in build_label_fields(species).items():
        data.setdefault(key, value)

    if len(parsed) == 1:
        part = parsed[0]
        data.setdefault("species_label", species[0])
        data.setdefault("species_formula", str(part["formula"]))
        data.setdefault("species_smiles", str(part["smiles"]))
        data.setdefault("species_wl", str(part["wl"]))
        data.setdefault("species_is_ads", bool(part["is_ads"]))
        data.setdefault("species_ads_n", int(part["ads_n"]))

    return node_id


def ensure_reaction_node(
    graph: nx.DiGraph,
    reactants_sorted: list[str],
    products_sorted: list[str],
    ttype: str,
) -> str:
    signature = reaction_signature(reactants_sorted, products_sorted, ttype)
    reaction_id = get_or_create_reaction_node_id(signature)
    if not graph.has_node(reaction_id):
        graph.add_node(reaction_id)

    data = graph.nodes[reaction_id]
    reactant_text = " + ".join(reactants_sorted)
    product_text = " + ".join(products_sorted)
    label = f"{reactant_text}\n->\n{product_text}"

    data.setdefault("orig_id", f"{reactant_text} -> {product_text}|type={ttype}")
    data.setdefault("label", label)
    data.setdefault("display_label", str(ttype or "reaction"))
    data.setdefault("node_type", "reaction")
    data.setdefault("is_reaction", True)
    data.setdefault("bipartite", "reaction")
    data.setdefault("reactants", reactant_text)
    data.setdefault("products", product_text)
    data.setdefault("type", str(ttype or ""))
    data.setdefault("reaction_reactants", json.dumps(reactants_sorted, ensure_ascii=False))
    data.setdefault("reaction_products", json.dumps(products_sorted, ensure_ascii=False))
    data.setdefault("reaction_equation", f"{reactant_text} -> {product_text}")
    data.setdefault("reaction_type", str(ttype or ""))
    data.setdefault("weight", 0)
    data.setdefault("n_reactants", len(reactants_sorted))
    data.setdefault("n_products", len(products_sorted))
    return reaction_id


def add_transform(
    graph: nx.DiGraph,
    reactants: list[str],
    products: list[str],
    ttype: str,
    frame: int | None = None,
) -> None:
    if not reactants or not products:
        return

    reactants_sorted = sorted(reactants)
    products_sorted = sorted(products)
    src = ensure_species_node(graph, " + ".join(reactants_sorted), reactants_sorted)
    dst = ensure_species_node(graph, " + ".join(products_sorted), products_sorted)

    if graph.has_edge(src, dst):
        graph[src][dst]["weight"] = int(graph[src][dst].get("weight", 0)) + 1
        old_type = graph[src][dst].get("type")
        if old_type is None:
            graph[src][dst]["type"] = ttype
        elif old_type != ttype:
            types = set(str(old_type).split("|"))
            types.add(ttype)
            graph[src][dst]["type"] = "|".join(sorted(types))
    else:
        graph.add_edge(src, dst, weight=1, type=ttype)

    edge = graph[src][dst]
    if frame is not None:
        frames = parse_frames_attr(edge.get("frames"))
        frames.append(int(frame))
        edge["frames"] = frames
    edge["src_is_ads"] = bool(graph.nodes[src].get("is_ads", False))
    edge["dst_is_ads"] = bool(graph.nodes[dst].get("is_ads", False))
    edge["is_ads_change"] = bool(edge["src_is_ads"] != edge["dst_is_ads"])
    edge["label"] = f"{edge.get('type')} ({edge.get('weight', 1)})"
    edge.setdefault("display_label", str(edge.get("type", "")))


def add_transform_bipartite(
    graph: nx.DiGraph,
    event_id: int,
    reactants: list[str],
    products: list[str],
    ttype: str,
) -> None:
    if not reactants or not products:
        return

    reactants_sorted = sorted(reactants)
    products_sorted = sorted(products)
    reaction_id = ensure_reaction_node(graph, reactants_sorted, products_sorted, ttype)
    graph.nodes[reaction_id]["weight"] = int(graph.nodes[reaction_id].get("weight", 0)) + 1

    for species, stoich in Counter(reactants_sorted).items():
        species_id = ensure_species_node(graph, species, [species])
        if graph.has_edge(species_id, reaction_id):
            graph[species_id][reaction_id]["weight"] = int(graph[species_id][reaction_id].get("weight", 0)) + 1
        else:
            graph.add_edge(species_id, reaction_id, role="reactant", stoich=int(stoich), weight=1)
        edge = graph[species_id][reaction_id]
        edge.setdefault("role", "reactant")
        edge.setdefault("stoich", int(stoich))
        edge["label"] = f"reactant x{edge.get('stoich', 1)} (n={edge.get('weight', 1)})"
        edge.setdefault("display_label", "reactant")

    for species, stoich in Counter(products_sorted).items():
        species_id = ensure_species_node(graph, species, [species])
        if graph.has_edge(reaction_id, species_id):
            graph[reaction_id][species_id]["weight"] = int(graph[reaction_id][species_id].get("weight", 0)) + 1
        else:
            graph.add_edge(reaction_id, species_id, role="product", stoich=int(stoich), weight=1)
        edge = graph[reaction_id][species_id]
        edge.setdefault("role", "product")
        edge.setdefault("stoich", int(stoich))
        edge["label"] = f"product x{edge.get('stoich', 1)} (n={edge.get('weight', 1)})"
        edge.setdefault("display_label", "product")


def summarize_reversible_reactions(
    graph: nx.DiGraph,
    use: str = "orig_id",
    min_total_weight: int = 1,
    sort_by: str = "total",
    include_frames: bool = False,
    unique_frames: bool = False,
) -> list[str]:
    def node_name(node_id) -> str:
        if use == "id":
            return str(node_id)
        return str(graph.nodes[node_id].get(use) or graph.nodes[node_id].get("label") or graph.nodes[node_id].get("orig_id") or node_id)

    edge_weights: dict[tuple[str, str], int] = {}
    edge_frames: dict[tuple[str, str], list[int]] = {}
    for u, v, data in graph.edges(data=True):
        try:
            weight = int(data.get("weight", 1))
        except Exception:
            try:
                weight = int(float(data.get("weight", 1)))
            except Exception:
                weight = 1
        edge_weights[(u, v)] = weight
        edge_frames[(u, v)] = parse_frames_attr(data.get("frames"))

    seen: set[frozenset] = set()
    rows: list[dict[str, object]] = []
    for u, v in graph.edges():
        pair = frozenset((u, v))
        if pair in seen:
            continue
        seen.add(pair)

        nu, nv = node_name(u), node_name(v)
        if (nu, nv) <= (nv, nu):
            a, b = u, v
            na, nb = nu, nv
        else:
            a, b = v, u
            na, nb = nv, nu

        w_fwd = edge_weights.get((a, b), 0)
        w_rev = edge_weights.get((b, a), 0)
        w_total = w_fwd + w_rev
        if w_total < min_total_weight:
            continue

        fwd_frames = list(edge_frames.get((a, b), []))
        rev_frames = list(edge_frames.get((b, a), []))
        if unique_frames:
            fwd_frames = sorted(set(fwd_frames))
            rev_frames = sorted(set(rev_frames))

        rows.append(
            {
                "na": na,
                "nb": nb,
                "w_fwd": w_fwd,
                "w_rev": w_rev,
                "w_total": w_total,
                "fwd_frames": fwd_frames,
                "rev_frames": rev_frames,
            }
        )

    if sort_by == "total":
        rows.sort(key=lambda row: int(row["w_total"]), reverse=True)
    elif sort_by == "name":
        rows.sort(key=lambda row: (str(row["na"]), str(row["nb"])))

    if include_frames:
        lines = ["reaction\tfwd\trev\ttotal\tfwd_frames\trev_frames"]
        for row in rows:
            lines.append(
                f"{row['na']}<->{row['nb']}\t{row['w_fwd']}\t{row['w_rev']}\t{row['w_total']}\t"
                f"{json.dumps(row['fwd_frames'], ensure_ascii=False)}\t"
                f"{json.dumps(row['rev_frames'], ensure_ascii=False)}"
            )
        return lines

    lines = ["reaction\tfwd\trev\ttotal"]
    for row in rows:
        lines.append(f"{row['na']}<->{row['nb']}\t{row['w_fwd']}\t{row['w_rev']}\t{row['w_total']}")
    return lines
