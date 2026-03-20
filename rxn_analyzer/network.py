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

_NODE_ID_MAP: dict[str, str] = {}
_NEXT_NODE_ID: int = 0
_RXN_NODE_ID_MAP: dict[tuple, str] = {}
_NEXT_RXN_ID: int = 0


def ensure_graph() -> nx.DiGraph:
    return nx.DiGraph()


def ensure_bipartite_graph() -> nx.DiGraph:
    return nx.DiGraph()


def reset_node_id_mapping() -> None:
    global _NODE_ID_MAP, _NEXT_NODE_ID, _RXN_NODE_ID_MAP, _NEXT_RXN_ID
    _NODE_ID_MAP = {}
    _NEXT_NODE_ID = 0
    _RXN_NODE_ID_MAP = {}
    _NEXT_RXN_ID = 0


def _parse_formula_counts(formula: str) -> Counter[str]:
    c: Counter[str] = Counter()
    for el, num in _ELEM_RE.findall(formula):
        c[el] += int(num) if num else 1
    return c


def _parse_species_label(s: str) -> dict[str, object]:
    formula = s.split("|", 1)[0]
    m = _SMILES_RE.search(s)
    smiles = m.group(1) if m else ""
    m = _WL_RE.search(s)
    wl = m.group(1) if m else ""
    is_ads = "|ads" in s
    m = _ADS_N_RE.search(s)
    ads_n = int(m.group(1)) if m else 0
    short = f"{formula}|ads" if is_ads else formula
    return {
        "formula": formula,
        "smiles": smiles,
        "wl": wl,
        "is_ads": is_ads,
        "ads_n": ads_n,
        "short": short,
    }


def _node_attrs_from_species_list(species: list[str]) -> dict[str, object]:
    formula_parts = [s.split("|", 1)[0] for s in species]
    elem_counts: Counter[str] = Counter()
    for f in formula_parts:
        elem_counts += _parse_formula_counts(f)

    ads_ns: list[int] = []
    is_ads_any = False
    for s in species:
        if "|ads" in s:
            is_ads_any = True
        m = _ADS_N_RE.search(s)
        if m:
            ads_ns.append(int(m.group(1)))

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


def _build_label_fields(species: list[str]) -> dict[str, str]:
    parts = [_parse_species_label(s) for s in species]
    label_formula = " + ".join(str(p["formula"]) for p in parts)
    label_smiles = " + ".join(str(p["smiles"]) for p in parts if p["smiles"])
    label_wl = " + ".join(str(p["wl"]) for p in parts if p["wl"])
    label_ads = "ads" if any(bool(p["is_ads"]) for p in parts) else "non_ads"
    label_ads_n = " + ".join(str(p["ads_n"]) for p in parts if bool(p["is_ads"]))
    label_short = " + ".join(str(p["short"]) for p in parts)
    return {
        "label_formula": label_formula,
        "label_smiles": label_smiles,
        "label_wl": label_wl,
        "label_ads": label_ads,
        "label_ads_n": label_ads_n,
        "label_short": label_short,
    }


def _get_or_create_numeric_node_id(long_id: str) -> str:
    global _NEXT_NODE_ID
    if long_id in _NODE_ID_MAP:
        return _NODE_ID_MAP[long_id]
    num_id = str(_NEXT_NODE_ID)
    _NODE_ID_MAP[long_id] = num_id
    _NEXT_NODE_ID += 1
    return num_id


def _ensure_node(G: nx.DiGraph, long_id: str, species: list[str]) -> str:
    nid = _get_or_create_numeric_node_id(long_id)
    if not G.has_node(nid):
        G.add_node(nid)

    nd = G.nodes[nid]
    parsed = [_parse_species_label(s) for s in species]
    short_display = " + ".join(str(p["short"]) for p in parsed) if parsed else long_id

    nd.setdefault("orig_id", long_id)
    nd.setdefault("label", long_id.replace(" + ", " +\n").replace("|ads", "\n|ads"))
    nd.setdefault("display_label", short_display)
    nd.setdefault("node_type", "species")
    nd.setdefault("is_reaction", False)
    nd.setdefault("bipartite", "species")
    nd.setdefault("species_labels", json.dumps(species, ensure_ascii=False))

    for k, v in _node_attrs_from_species_list(species).items():
        nd.setdefault(k, v)
    for k, v in _build_label_fields(species).items():
        nd.setdefault(k, v)

    if len(parsed) == 1:
        p = parsed[0]
        nd.setdefault("species_label", species[0])
        nd.setdefault("species_formula", str(p["formula"]))
        nd.setdefault("species_smiles", str(p["smiles"]))
        nd.setdefault("species_wl", str(p["wl"]))
        nd.setdefault("species_is_ads", bool(p["is_ads"]))
        nd.setdefault("species_ads_n", int(p["ads_n"]))

    return nid


def _reaction_signature(reactants_sorted: list[str], products_sorted: list[str], ttype: str) -> tuple:
    return (tuple(reactants_sorted), tuple(products_sorted), str(ttype or ""))


def _get_or_create_reaction_node_id(sig: tuple) -> str:
    global _NEXT_RXN_ID
    if sig in _RXN_NODE_ID_MAP:
        return _RXN_NODE_ID_MAP[sig]
    rid = f"r{_NEXT_RXN_ID}"
    _RXN_NODE_ID_MAP[sig] = rid
    _NEXT_RXN_ID += 1
    return rid


def _ensure_reaction_node(
    G: nx.DiGraph,
    reactants_sorted: list[str],
    products_sorted: list[str],
    ttype: str,
) -> str:
    sig = _reaction_signature(reactants_sorted, products_sorted, ttype)
    rid = _get_or_create_reaction_node_id(sig)
    if not G.has_node(rid):
        G.add_node(rid)

    nd = G.nodes[rid]
    r_str = " + ".join(reactants_sorted)
    p_str = " + ".join(products_sorted)
    label = f"{r_str}\n->\n{p_str}"

    nd.setdefault("orig_id", f"{r_str} -> {p_str}|type={ttype}")
    nd.setdefault("label", label)
    nd.setdefault("display_label", str(ttype or "reaction"))
    nd.setdefault("node_type", "reaction")
    nd.setdefault("is_reaction", True)
    nd.setdefault("bipartite", "reaction")
    nd.setdefault("reactants", r_str)
    nd.setdefault("products", p_str)
    nd.setdefault("type", str(ttype or ""))
    nd.setdefault("reaction_reactants", json.dumps(reactants_sorted, ensure_ascii=False))
    nd.setdefault("reaction_products", json.dumps(products_sorted, ensure_ascii=False))
    nd.setdefault("reaction_equation", f"{r_str} -> {p_str}")
    nd.setdefault("reaction_type", str(ttype or ""))
    nd.setdefault("weight", 0)
    nd.setdefault("n_reactants", len(reactants_sorted))
    nd.setdefault("n_products", len(products_sorted))

    return rid


def _parse_frames_attr(v) -> list[int]:
    if v is None or isinstance(v, bool):
        return []
    if isinstance(v, (int, float)):
        try:
            return [int(v)]
        except Exception:
            return []
    if isinstance(v, (list, tuple, set)):
        out: list[int] = []
        for x in v:
            try:
                out.append(int(x))
            except Exception:
                pass
        return out
    if isinstance(v, str):
        s = v.strip()
        if not s:
            return []
        for parser in (json.loads, ast.literal_eval):
            try:
                obj = parser(s)
            except Exception:
                continue
            if isinstance(obj, (list, tuple, set)):
                out: list[int] = []
                for x in obj:
                    try:
                        out.append(int(x))
                    except Exception:
                        pass
                return out
            try:
                return [int(obj)]
            except Exception:
                pass
        s2 = s.strip("[](){}")
        if not s2:
            return []
        out = []
        for tok in re.split(r"[,\s]+", s2):
            if not tok:
                continue
            try:
                out.append(int(tok))
            except Exception:
                pass
        return out
    return []


def add_transform(
    G: nx.DiGraph,
    reactants: list[str],
    products: list[str],
    ttype: str,
    frame: int | None = None,
):
    if not reactants or not products:
        return

    reactants_sorted = sorted(reactants)
    products_sorted = sorted(products)
    long_src = " + ".join(reactants_sorted)
    long_dst = " + ".join(products_sorted)

    src = _ensure_node(G, long_src, reactants_sorted)
    dst = _ensure_node(G, long_dst, products_sorted)

    if G.has_edge(src, dst):
        G[src][dst]["weight"] = int(G[src][dst].get("weight", 0)) + 1
        old = G[src][dst].get("type")
        if old is None:
            G[src][dst]["type"] = ttype
        elif old != ttype:
            types = set(str(old).split("|"))
            types.add(ttype)
            G[src][dst]["type"] = "|".join(sorted(types))
    else:
        G.add_edge(src, dst, weight=1, type=ttype)

    ed = G[src][dst]
    if frame is not None:
        f_list = _parse_frames_attr(ed.get("frames"))
        f_list.append(int(frame))
        ed["frames"] = f_list
    ed["src_is_ads"] = bool(G.nodes[src].get("is_ads", False))
    ed["dst_is_ads"] = bool(G.nodes[dst].get("is_ads", False))
    ed["is_ads_change"] = bool(ed["src_is_ads"] != ed["dst_is_ads"])
    ed["label"] = f"{ed.get('type')} ({ed.get('weight', 1)})"
    ed.setdefault("display_label", str(ed.get("type", "")))


def add_transform_bipartite(
    G: nx.DiGraph,
    event_id: int,
    reactants: list[str],
    products: list[str],
    ttype: str,
):
    if not reactants or not products:
        return

    reactants_sorted = sorted(reactants)
    products_sorted = sorted(products)
    rid = _ensure_reaction_node(G, reactants_sorted, products_sorted, ttype)
    G.nodes[rid]["weight"] = int(G.nodes[rid].get("weight", 0)) + 1

    r_counts = Counter(reactants_sorted)
    for sp, stoich in r_counts.items():
        sid = _ensure_node(G, sp, [sp])
        if G.has_edge(sid, rid):
            G[sid][rid]["weight"] = int(G[sid][rid].get("weight", 0)) + 1
        else:
            G.add_edge(sid, rid, role="reactant", stoich=int(stoich), weight=1)
        ed = G[sid][rid]
        ed.setdefault("role", "reactant")
        ed.setdefault("stoich", int(stoich))
        ed["label"] = f"reactant x{ed.get('stoich', 1)} (n={ed.get('weight', 1)})"
        ed.setdefault("display_label", "reactant")

    p_counts = Counter(products_sorted)
    for sp, stoich in p_counts.items():
        sid = _ensure_node(G, sp, [sp])
        if G.has_edge(rid, sid):
            G[rid][sid]["weight"] = int(G[rid][sid].get("weight", 0)) + 1
        else:
            G.add_edge(rid, sid, role="product", stoich=int(stoich), weight=1)
        ed = G[rid][sid]
        ed.setdefault("role", "product")
        ed.setdefault("stoich", int(stoich))
        ed["label"] = f"product x{ed.get('stoich', 1)} (n={ed.get('weight', 1)})"
        ed.setdefault("display_label", "product")


def finalize_graph_for_export(G: nx.DiGraph) -> nx.DiGraph:
    def scalarize(v):
        if v is None:
            return ""
        if isinstance(v, (str, int, float, bool)):
            return v
        return str(v)

    for k in list(G.graph.keys()):
        G.graph[k] = scalarize(G.graph[k])
    for _n, d in G.nodes(data=True):
        for k in list(d.keys()):
            d[k] = scalarize(d[k])
    for _u, _v, d in G.edges(data=True):
        for k in list(d.keys()):
            d[k] = scalarize(d[k])
    return G


def summarize_reversible_reactions(
    G: nx.DiGraph,
    use: str = "orig_id",
    min_total_weight: int = 1,
    sort_by: str = "total",
    include_frames: bool = False,
    unique_frames: bool = False,
) -> list[str]:
    def node_name(n) -> str:
        if use == "id":
            return str(n)
        return str(G.nodes[n].get(use) or G.nodes[n].get("label") or G.nodes[n].get("orig_id") or n)

    w_dir: dict[tuple[str, str], int] = {}
    frames_dir: dict[tuple[str, str], list[int]] = {}
    for u, v, d in G.edges(data=True):
        try:
            w = int(d.get("weight", 1))
        except Exception:
            try:
                w = int(float(d.get("weight", 1)))
            except Exception:
                w = 1
        w_dir[(u, v)] = w
        frames_dir[(u, v)] = _parse_frames_attr(d.get("frames"))

    seen: set[frozenset] = set()
    rows: list[dict[str, object]] = []
    for u, v in G.edges():
        key = frozenset((u, v))
        if key in seen:
            continue
        seen.add(key)

        nu, nv = node_name(u), node_name(v)
        if (nu, nv) <= (nv, nu):
            a, b = u, v
            na, nb = nu, nv
        else:
            a, b = v, u
            na, nb = nv, nu

        w_fwd = w_dir.get((a, b), 0)
        w_rev = w_dir.get((b, a), 0)
        w_total = w_fwd + w_rev
        if w_total < min_total_weight:
            continue

        fwd_frames = list(frames_dir.get((a, b), []))
        rev_frames = list(frames_dir.get((b, a), []))
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
        rows.sort(key=lambda x: int(x["w_total"]), reverse=True)
    elif sort_by == "name":
        rows.sort(key=lambda x: (str(x["na"]), str(x["nb"])))

    if include_frames:
        lines = ["reaction\tfwd\trev\ttotal\tfwd_frames\trev_frames"]
        for r in rows:
            lines.append(
                f"{r['na']}<->{r['nb']}\t{r['w_fwd']}\t{r['w_rev']}\t{r['w_total']}\t"
                f"{json.dumps(r['fwd_frames'], ensure_ascii=False)}\t"
                f"{json.dumps(r['rev_frames'], ensure_ascii=False)}"
            )
    else:
        lines = ["reaction\tfwd\trev\ttotal"]
        for r in rows:
            lines.append(f"{r['na']}<->{r['nb']}\t{r['w_fwd']}\t{r['w_rev']}\t{r['w_total']}")
    return lines
