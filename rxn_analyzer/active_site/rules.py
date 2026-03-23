from __future__ import annotations

from collections import Counter

from ase import Atoms


def formula_from_atoms(atoms: Atoms, indices: list[int] | tuple[int, ...]) -> str:
    if not indices:
        return "bare"

    if len(indices) == 2:
        syms2 = sorted(atoms[i].symbol for i in indices)
        if syms2 == ["H", "O"]:
            return "OH"

    syms = [atoms[i].symbol for i in indices]
    counts = Counter(syms)

    keys: list[str] = []
    if "C" in counts:
        keys.append("C")
    if "H" in counts:
        keys.append("H")
    for key in sorted(counts):
        if key not in {"C", "H"}:
            keys.append(key)

    out: list[str] = []
    for key in keys:
        n = counts[key]
        out.append(key if n == 1 else f"{key}{n}")
    return "".join(out) if out else "bare"


def topology_from_state_components(
    atoms: Atoms,
    component_nodes: list[list[int]],
    *,
    core_contacts: set[int],
) -> str:
    if not component_nodes:
        return "bare"

    tokens: list[str] = []
    for comp in component_nodes:
        f = formula_from_atoms(atoms, comp)
        if f == "HO" and core_contacts.intersection(comp):
            tokens.append("OH")
        elif f == "H2O" and core_contacts.intersection(comp):
            tokens.append("H2O")
        else:
            tokens.append(f)

    counts = Counter(tokens)
    ordered = sorted(counts)
    out: list[str] = []
    for token in ordered:
        n = counts[token]
        if token == "OH":
            out.append("(OH)" if n == 1 else f"(OH){n}")
        else:
            out.append(token if n == 1 else f"{token}{n}")
    return "-".join(out) if out else "bare"


def _build_topo_components(member_list: list[int], cov_adj: dict[int, set[int]]) -> list[list[int]]:
    topo_components: list[list[int]] = []
    visited: set[int] = set()
    member_set = set(member_list)
    for atom in member_list:
        if atom in visited:
            continue
        stack = [atom]
        visited.add(atom)
        comp = [atom]
        while stack:
            x = stack.pop()
            for y in cov_adj.get(x, set()):
                if y in member_set and y not in visited:
                    visited.add(y)
                    stack.append(y)
                    comp.append(y)
        topo_components.append(sorted(comp))
    return topo_components


def _pick_primary_core_metal(
    atoms: Atoms,
    core_members: list[int] | tuple[int, ...],
    *,
    excluded: set[str] | None = None,
) -> tuple[str | None, set[int]]:
    excluded = excluded or set()
    syms = atoms.get_chemical_symbols()
    metals = {
        atom
        for atom in core_members
        if syms[atom] not in excluded and syms[atom] not in {"H", "C", "N", "O", "Si", "P", "S"}
    }
    if not metals:
        return None, set()
    ordered = sorted(metals)
    primary = syms[ordered[0]].lower()
    same_symbol = {atom for atom in ordered if syms[atom].lower() == primary}
    return primary, same_symbol


def _format_metal_base_label(
    metal_label: str,
    *,
    n_mh: int,
    n_moh: int,
    n_mo: int,
    has_carbon: bool,
) -> str:
    parts: list[str] = []
    if n_mh > 0:
        parts.append("h" if n_mh == 1 else f"h{n_mh}")
    if n_moh > 0:
        parts.append("oh" if n_moh == 1 else f"oh{n_moh}")
    if n_mo > 0:
        parts.append("o" if n_mo == 1 else f"o{n_mo}")
    if has_carbon:
        parts.append("organic")
    if not parts:
        return f"{metal_label}_bare"
    return f"{metal_label}_" + "_".join(parts)


def _classify_bas_state(
    *,
    atoms: Atoms,
    site_family: str,
    core_members: list[int] | tuple[int, ...],
    state_body: list[int],
    associated_topology: str,
    state_formula: str,
    state_topology: str,
    cov_adj: dict[int, set[int]],
    core_edge_map: dict[int, set[int]],
) -> tuple[str, dict[str, object]]:
    syms = atoms.get_chemical_symbols()
    oh_count = 0
    direct_h_count = 0
    direct_o_count = 0

    for atom in state_body:
        core_neighbors = set(core_edge_map.get(atom, set())) & set(core_members)
        if not core_neighbors:
            continue
        if syms[atom] == "H":
            direct_h_count += 1
        elif syms[atom] == "O":
            direct_o_count += 1
            h_neighbors = {nbr for nbr in cov_adj.get(atom, set()) if nbr in state_body and syms[nbr] == "H"}
            if h_neighbors:
                oh_count += 1

    has_incorporated_carbon = any(syms[atom] == "C" for atom in state_body)
    has_associated_water = "H2O" in associated_topology
    family_lc = site_family.lower()

    if has_incorporated_carbon:
        base_label = "alkoxy_bas"
    elif direct_o_count > 0 and direct_h_count == 0 and oh_count == 0:
        base_label = "deprotonated_bas"
    elif oh_count > 0 or direct_h_count > 0:
        base_label = "protonated_bas"
    elif state_topology in {"", "bare"} and state_formula == "bare":
        base_label = "bare_bas"
    else:
        base_label = "unknown_bas_state"

    if "ga" in family_lc:
        base_label = base_label.replace("_bas", "_bas_ga")

    descriptors = {
        "n_bas_oh": oh_count,
        "n_bas_direct_h": direct_h_count,
        "n_bas_direct_o": direct_o_count,
        "has_incorporated_carbon": has_incorporated_carbon,
        "has_associated_water": has_associated_water,
    }

    if has_associated_water:
        return f"{base_label} + assoc(H2O)", descriptors
    if associated_topology not in {"", "bare"}:
        return f"{base_label} + assoc({associated_topology})", descriptors
    return base_label, descriptors


def classify_active_site_state(
    *,
    atoms: Atoms,
    site_family: str,
    core_members: list[int] | tuple[int, ...],
    intrinsic_members: list[int] | tuple[int, ...],
    incorporated_members: list[int] | tuple[int, ...],
    associated_members: list[int] | tuple[int, ...],
    associated_topology: str,
    state_formula: str,
    state_topology: str,
    cov_adj: dict[int, set[int]],
    core_edge_map: dict[int, set[int]],
) -> tuple[str, dict[str, object]]:
    family_lc = site_family.lower()
    state_body = sorted(set(intrinsic_members) | set(incorporated_members))
    state_syms = atoms.get_chemical_symbols()

    if "las" in family_lc:
        metal_label, metal_core = _pick_primary_core_metal(atoms, core_members)
        if metal_label and metal_core:
            direct_h = {
                atom
                for atom in state_body
                if state_syms[atom] == "H" and bool(set(core_edge_map.get(atom, set())) & metal_core)
            }
            direct_o = {
                atom
                for atom in state_body
                if state_syms[atom] == "O" and bool(set(core_edge_map.get(atom, set())) & metal_core)
            }

            n_moh = 0
            n_mo = 0
            for o_atom in sorted(direct_o):
                h_neighbors = {
                    nbr
                    for nbr in cov_adj.get(o_atom, set())
                    if nbr in state_body and state_syms[nbr] == "H"
                }
                if h_neighbors:
                    n_moh += 1
                else:
                    n_mo += 1

            # Hydrides are H directly bound to the metal and not already counted as hydroxyl H.
            metal_bound_oh_h = {
                nbr
                for o_atom in direct_o
                for nbr in cov_adj.get(o_atom, set())
                if nbr in state_body and state_syms[nbr] == "H"
            }
            n_mh = len(direct_h - metal_bound_oh_h)
            has_carbon = any(state_syms[atom] == "C" for atom in state_body)

            base_label = _format_metal_base_label(
                metal_label,
                n_mh=n_mh,
                n_moh=n_moh,
                n_mo=n_mo,
                has_carbon=has_carbon,
            )

            descriptors = {
                "metal_center_symbol": metal_label.upper(),
                "n_metal_h": n_mh,
                "n_metal_oh": n_moh,
                "n_metal_o": n_mo,
                "has_incorporated_carbon": has_carbon,
            }
            descriptors[f"n_{metal_label}_mh"] = n_mh
            descriptors[f"n_{metal_label}_moh"] = n_moh
            descriptors[f"n_{metal_label}_mo"] = n_mo

            if associated_topology not in {"", "bare"}:
                return f"{base_label} + assoc({associated_topology})", descriptors
            return base_label, descriptors

    if "bas" in family_lc:
        return _classify_bas_state(
            atoms=atoms,
            site_family=site_family,
            core_members=core_members,
            state_body=state_body,
            associated_topology=associated_topology,
            state_formula=state_formula,
            state_topology=state_topology,
            cov_adj=cov_adj,
            core_edge_map=core_edge_map,
        )

    base_label = state_topology if state_topology not in {"", "bare"} else state_formula
    if associated_topology not in {"", "bare"}:
        return f"{base_label if base_label not in {'', 'bare'} else 'bare'} + assoc({associated_topology})", {}
    return base_label if base_label not in {"", "bare"} else "bare", {}
