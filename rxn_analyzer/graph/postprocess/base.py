"""Shared helpers for graph postprocessing."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from typing import Any

import networkx as nx


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(value)
    except Exception:
        return default


def summarize_graph(graph: nx.Graph, tag: str) -> None:
    print(f"[{tag}] nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")


def edge_weight(data: dict[str, Any]) -> int:
    return safe_int(data.get("weight", 1), 1)


def filter_edges_by_min_weight(graph: nx.Graph, min_weight: int) -> nx.Graph:
    if min_weight <= 1:
        return graph.copy()
    out = graph.copy()
    for u, v, data in list(out.edges(data=True)):
        if edge_weight(data) < min_weight:
            out.remove_edge(u, v)
    return out


def prune_isolates(graph: nx.Graph) -> nx.Graph:
    out = graph.copy()
    out.remove_nodes_from(list(nx.isolates(out)))
    return out


def _is_reaction_bipartite_graph(graph: nx.Graph) -> bool:
    node_types = {str(data.get("node_type", "")) for _node_id, data in graph.nodes(data=True)}
    return "reaction" in node_types and "species" in node_types


def _parse_species_list_attr(value: object) -> tuple[str, ...]:
    if value is None:
        return tuple()
    if isinstance(value, (list, tuple)):
        return tuple(str(x) for x in value)
    text = str(value).strip()
    if not text:
        return tuple()
    try:
        parsed = json.loads(text)
        if isinstance(parsed, list):
            return tuple(str(x) for x in parsed)
    except Exception:
        pass
    return tuple(x.strip() for x in text.split(" + ") if x.strip())


def _collapse_reversible_bipartite(graph: nx.Graph) -> nx.Graph:
    out = nx.DiGraph()
    species_node_by_label: dict[str, str] = {}

    for node_id, data in graph.nodes(data=True):
        if str(data.get("node_type", "")) == "species":
            out.add_node(node_id, **data)
            species_node_by_label[str(data.get("orig_id", ""))] = str(node_id)

    active_reaction_nodes = {
        str(node_id)
        for node_id, data in graph.nodes(data=True)
        if str(data.get("node_type", "")) == "reaction" and graph.degree(node_id) > 0
    }

    grouped: dict[tuple[tuple[str, ...], tuple[str, ...]], list[tuple[str, dict[str, Any], tuple[str, ...], tuple[str, ...]]]] = defaultdict(list)
    for node_id, data in graph.nodes(data=True):
        if str(data.get("node_type", "")) != "reaction":
            continue
        if str(node_id) not in active_reaction_nodes:
            continue
        reactants = tuple(sorted(_parse_species_list_attr(data.get("reaction_reactants"))))
        products = tuple(sorted(_parse_species_list_attr(data.get("reaction_products"))))
        canonical = tuple(sorted((reactants, products)))
        grouped[canonical].append((str(node_id), data, reactants, products))

    next_index = 0

    def add_or_merge_edge(
        source: str,
        target: str,
        *,
        role: str,
        side: str,
        stoich: int,
        weight: int,
        display_label: str,
    ) -> None:
        if out.has_edge(source, target):
            out[source][target]["weight"] = safe_int(out[source][target].get("weight", 0), 0) + weight
            out[source][target]["stoich"] = max(safe_int(out[source][target].get("stoich", 0), 0), stoich)
        else:
            out.add_edge(
                source,
                target,
                role=role,
                side=side,
                stoich=stoich,
                weight=weight,
                label=f"{display_label} (n={weight})",
                display_label=display_label,
            )

    for canonical, members in grouped.items():
        if not members:
            continue

        if len(canonical) == 2:
            side_a, side_b = canonical
        else:
            side_a, side_b = tuple(), tuple()

        direction_set = {(reactants, products) for _node_id, _data, reactants, products in members}
        has_ab = (side_a, side_b) in direction_set
        has_ba = (side_b, side_a) in direction_set
        reversible = bool(has_ab and has_ba and side_a != side_b)

        if reversible:
            lhs, rhs = side_a, side_b
        else:
            # Keep the dominant observed direction for one-way reactions.
            direction_weight: dict[tuple[tuple[str, ...], tuple[str, ...]], int] = defaultdict(int)
            for _old_node_id, data, reactants, products in members:
                direction_weight[(reactants, products)] += safe_int(data.get("weight", 1), 1)
            lhs, rhs = max(
                direction_weight.items(),
                key=lambda item: (item[1], item[0][0], item[0][1]),
            )[0]

        merged_id = f"rr{next_index}"
        next_index += 1

        node_weight = 0
        forward_weight = 0
        reverse_weight = 0
        type_tokens: set[str] = set()
        source_nodes: list[str] = []

        for old_node_id, data, reactants, products in members:
            source_nodes.append(old_node_id)
            node_weight += safe_int(data.get("weight", 1), 1)
            type_tokens.update(token.strip() for token in str(data.get("reaction_type", data.get("type", ""))).split("|") if token.strip())
            if reactants == lhs and products == rhs:
                forward_weight += safe_int(data.get("weight", 1), 1)
            elif reactants == rhs and products == lhs:
                reverse_weight += safe_int(data.get("weight", 1), 1)
            elif not reversible:
                forward_weight += safe_int(data.get("weight", 1), 1)

        equation = (
            f"{' + '.join(lhs)} <-> {' + '.join(rhs)}"
            if reversible
            else f"{' + '.join(lhs)} -> {' + '.join(rhs)}"
        )
        display = "<->" if reversible else "reaction"
        out.add_node(
            merged_id,
            node_type="reaction",
            is_reaction=True,
            bipartite="reaction",
            label=equation.replace(" <-> ", "\n<->\n").replace(" -> ", "\n->\n"),
            display_label=display,
            orig_id=equation,
            reaction_reactants=json.dumps(lhs, ensure_ascii=False),
            reaction_products=json.dumps(rhs, ensure_ascii=False),
            reaction_equation=equation,
            reaction_type="|".join(sorted(type_tokens)),
            type="|".join(sorted(type_tokens)),
            weight=node_weight,
            weight_forward=forward_weight,
            weight_reverse=reverse_weight,
            weight_total=node_weight,
            is_reversible=reversible,
            source_reaction_nodes=json.dumps(source_nodes, ensure_ascii=False),
        )

        participant_species = sorted(set(lhs) | set(rhs))
        for species_label in participant_species:
            species_node_id = species_node_by_label.get(species_label)
            if species_node_id is None:
                continue

            lhs_stoich = lhs.count(species_label)
            rhs_stoich = rhs.count(species_label)
            if lhs_stoich > 0:
                add_or_merge_edge(
                    species_node_id,
                    merged_id,
                    role="reactant",
                    side="lhs",
                    stoich=lhs_stoich,
                    weight=node_weight,
                    display_label="reactant",
                )
            if rhs_stoich > 0:
                add_or_merge_edge(
                    merged_id,
                    species_node_id,
                    role="product",
                    side="rhs",
                    stoich=rhs_stoich,
                    weight=node_weight,
                    display_label="product",
                )
            if reversible:
                if rhs_stoich > 0:
                    add_or_merge_edge(
                        species_node_id,
                        merged_id,
                        role="reactant_reverse",
                        side="rhs_reverse",
                        stoich=rhs_stoich,
                        weight=reverse_weight or node_weight,
                        display_label="reactant_rev",
                    )
                if lhs_stoich > 0:
                    add_or_merge_edge(
                        merged_id,
                        species_node_id,
                        role="product_reverse",
                        side="lhs_reverse",
                        stoich=lhs_stoich,
                        weight=reverse_weight or node_weight,
                        display_label="product_rev",
                    )

    return out


def collapse_reversible(graph: nx.Graph) -> nx.Graph:
    if not _is_reaction_bipartite_graph(graph):
        raise ValueError("Graph postprocess collapse_reversible now only supports bipartite species/reaction graphs.")
    return _collapse_reversible_bipartite(graph)


def node_display_name(node_id: str, data: dict[str, Any]) -> str:
    return str(
        data.get("display_label")
        or data.get("label")
        or data.get("orig_id")
        or node_id
    )


def node_name_to_id(graph: nx.Graph, name: str) -> str | None:
    if graph.has_node(name):
        return str(name)
    for node_id, data in graph.nodes(data=True):
        if str(data.get("orig_id", "")) == name:
            return str(node_id)
        if str(data.get("label", "")) == name:
            return str(node_id)
        if str(data.get("display_label", "")) == name:
            return str(node_id)
    return None


def weighted_degree_scores(graph: nx.Graph) -> dict[str, int]:
    scores: dict[str, int] = defaultdict(int)
    if graph.is_directed():
        for node_id in graph.nodes():
            score = 0
            for _u, _v, data in graph.in_edges(node_id, data=True):
                score += edge_weight(data)
            for _u, _v, data in graph.out_edges(node_id, data=True):
                score += edge_weight(data)
            scores[str(node_id)] = score
    else:
        for node_id in graph.nodes():
            score = 0
            for _u, _v, data in graph.edges(node_id, data=True):
                score += edge_weight(data)
            scores[str(node_id)] = score
    return scores


def reaction_weight_scores(graph: nx.Graph) -> dict[str, int]:
    scores: dict[str, int] = defaultdict(int)
    fallback_scores = weighted_degree_scores(graph)

    for node_id, data in graph.nodes(data=True):
        node_type = str(data.get("node_type", ""))
        if node_type == "reaction":
            scores[str(node_id)] = safe_int(data.get("weight_total", data.get("weight", 0)), 0)
        elif node_type == "species":
            score = 0
            if graph.is_directed():
                neighbors = set(graph.predecessors(node_id)) | set(graph.neighbors(node_id))
            else:
                neighbors = set(graph.neighbors(node_id))
            for neighbor in neighbors:
                neighbor_data = graph.nodes[neighbor]
                if str(neighbor_data.get("node_type", "")) != "reaction":
                    continue
                score += safe_int(neighbor_data.get("weight_total", neighbor_data.get("weight", 0)), 0)
            scores[str(node_id)] = score
        else:
            scores[str(node_id)] = fallback_scores.get(str(node_id), 0)
    return scores


def compute_focus_scores(graph: nx.Graph, mode: str) -> dict[str, int]:
    mode_norm = str(mode or "weighted_degree").strip().lower()
    if mode_norm == "weighted_degree":
        return weighted_degree_scores(graph)
    if mode_norm == "reaction_weight":
        return reaction_weight_scores(graph)
    if mode_norm == "hybrid":
        degree_scores = weighted_degree_scores(graph)
        reaction_scores = reaction_weight_scores(graph)
        return {
            node_id: int(degree_scores.get(node_id, 0)) + int(reaction_scores.get(node_id, 0))
            for node_id in {str(n) for n in graph.nodes()}
        }
    raise ValueError(
        f"Unknown focus score_mode: {mode}. Supported modes: weighted_degree, reaction_weight, hybrid."
    )


def write_focus_scores_csv(
    path: str,
    graph: nx.Graph,
    scores: dict[str, int],
    *,
    score_mode: str,
    selected: set[str] | None = None,
) -> None:
    print(f"[write] {path}")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "node_id",
                "node_type",
                "display_name",
                "orig_id",
                "score_mode",
                "importance_score",
                "selected",
            ],
        )
        writer.writeheader()
        rows: list[tuple[int, str, dict[str, Any]]] = []
        for node_id, data in graph.nodes(data=True):
            rows.append((scores.get(str(node_id), 0), str(node_id), data))
        rows.sort(key=lambda item: (-item[0], item[1]))
        for score, node_id, data in rows:
            writer.writerow(
                {
                    "node_id": node_id,
                    "node_type": data.get("node_type", ""),
                    "display_name": node_display_name(node_id, data),
                    "orig_id": data.get("orig_id", ""),
                    "score_mode": score_mode,
                    "importance_score": score,
                    "selected": bool(selected and node_id in selected),
                }
            )
