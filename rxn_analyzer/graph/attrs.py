from __future__ import annotations

import ast
import json
import re
from collections import Counter

import networkx as nx


_ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
_ADS_N_RE = re.compile(r"\|ads\(n=(\d+)")
_SMILES_RE = re.compile(r"\|smiles=([^|]+)")
_WL_RE = re.compile(r"\|wl=([^|]+)")


def parse_formula_counts(formula: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    for el, num in _ELEM_RE.findall(formula):
        counts[el] += int(num) if num else 1
    return counts


def parse_species_label(label: str) -> dict[str, object]:
    formula = label.split("|", 1)[0]
    smiles_match = _SMILES_RE.search(label)
    wl_match = _WL_RE.search(label)
    ads_match = _ADS_N_RE.search(label)
    smiles = smiles_match.group(1) if smiles_match else ""
    wl = wl_match.group(1) if wl_match else ""
    is_ads = "|ads" in label
    ads_n = int(ads_match.group(1)) if ads_match else 0
    short = f"{formula}|ads" if is_ads else formula
    return {
        "formula": formula,
        "smiles": smiles,
        "wl": wl,
        "is_ads": is_ads,
        "ads_n": ads_n,
        "short": short,
    }


def node_attrs_from_species_list(species: list[str]) -> dict[str, object]:
    formula_parts = [s.split("|", 1)[0] for s in species]
    elem_counts: Counter[str] = Counter()
    for formula in formula_parts:
        elem_counts += parse_formula_counts(formula)

    ads_ns: list[int] = []
    is_ads_any = False
    for label in species:
        if "|ads" in label:
            is_ads_any = True
        ads_match = _ADS_N_RE.search(label)
        if ads_match:
            ads_ns.append(int(ads_match.group(1)))

    attrs: dict[str, object] = {
        "n_components": len(species),
        "is_mixture": len(species) > 1,
        "is_ads": bool(is_ads_any),
        "ads_n_total": int(sum(ads_ns)) if ads_ns else 0,
        "ads_n_max": int(max(ads_ns)) if ads_ns else 0,
        "nC": int(elem_counts.get("C", 0)),
        "nH": int(elem_counts.get("H", 0)),
        "nO": int(elem_counts.get("O", 0)),
        "nN": int(elem_counts.get("N", 0)),
        "nS": int(elem_counts.get("S", 0)),
    }
    for el, n in elem_counts.items():
        key = f"n{el}"
        if key not in attrs:
            attrs[key] = int(n)
    return attrs


def build_label_fields(species: list[str]) -> dict[str, str]:
    parts = [parse_species_label(s) for s in species]
    label_formula = " + ".join(str(part["formula"]) for part in parts)
    label_smiles = " + ".join(str(part["smiles"]) for part in parts if part["smiles"])
    label_wl = " + ".join(str(part["wl"]) for part in parts if part["wl"])
    label_ads = "ads" if any(bool(part["is_ads"]) for part in parts) else "non_ads"
    label_ads_n = " + ".join(str(part["ads_n"]) for part in parts if bool(part["is_ads"]))
    label_short = " + ".join(str(part["short"]) for part in parts)
    return {
        "label_formula": label_formula,
        "label_smiles": label_smiles,
        "label_wl": label_wl,
        "label_ads": label_ads,
        "label_ads_n": label_ads_n,
        "label_short": label_short,
    }


def parse_frames_attr(value) -> list[int]:
    if value is None or isinstance(value, bool):
        return []
    if isinstance(value, (int, float)):
        try:
            return [int(value)]
        except Exception:
            return []
    if isinstance(value, (list, tuple, set)):
        out: list[int] = []
        for item in value:
            try:
                out.append(int(item))
            except Exception:
                pass
        return out
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        for parser in (json.loads, ast.literal_eval):
            try:
                obj = parser(text)
            except Exception:
                continue
            if isinstance(obj, (list, tuple, set)):
                out: list[int] = []
                for item in obj:
                    try:
                        out.append(int(item))
                    except Exception:
                        pass
                return out
            try:
                return [int(obj)]
            except Exception:
                pass
        text = text.strip("[](){}")
        if not text:
            return []
        out = []
        for token in re.split(r"[,\s]+", text):
            if not token:
                continue
            try:
                out.append(int(token))
            except Exception:
                pass
        return out
    return []


def finalize_graph_for_export(graph: nx.DiGraph) -> nx.DiGraph:
    def scalarize(value):
        if value is None:
            return ""
        if isinstance(value, (str, int, float, bool)):
            return value
        return str(value)

    for key in list(graph.graph.keys()):
        graph.graph[key] = scalarize(graph.graph[key])
    for _node, data in graph.nodes(data=True):
        for key in list(data.keys()):
            data[key] = scalarize(data[key])
    for _u, _v, data in graph.edges(data=True):
        for key in list(data.keys()):
            data[key] = scalarize(data[key])
    return graph
