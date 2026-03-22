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
