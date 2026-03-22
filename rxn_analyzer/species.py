from __future__ import annotations
from dataclasses import dataclass
from collections import Counter, defaultdict
from abc import ABC, abstractmethod
import numpy as np
from ase import Atoms

from .edges import EdgeType


def connected_components(n_nodes: int, edges: set[tuple[int, int]], allowed: set[int]) -> list[list[int]]:
    adj = {u: [] for u in allowed}
    for i, j in edges:
        if i in allowed and j in allowed:
            adj[i].append(j)
            adj[j].append(i)
    seen = set()
    comps = []
    for u in sorted(allowed):
        if u in seen:
            continue
        stack = [u]
        seen.add(u)
        comp = [u]
        while stack:
            x = stack.pop()
            for y in adj.get(x, []):
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
                    comp.append(y)
        comps.append(sorted(comp))
    return comps


def wl_hash(atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]], iters: int = 3) -> str:
    node_set = set(nodes)
    adj = {u: [] for u in nodes}
    for i, j in cov_edges:
        if i in node_set and j in node_set:
            adj[i].append(j)
            adj[j].append(i)

    z = atoms.numbers
    labels = {u: f"Z{int(z[u])}" for u in nodes}

    for _ in range(iters):
        new_labels = {}
        for u in nodes:
            neigh = sorted(labels[v] for v in adj[u])
            s = labels[u] + "|" + ",".join(neigh)
            h = 1469598103934665603
            for ch in s.encode("utf-8"):
                h ^= ch
                h = (h * 1099511628211) & ((1 << 64) - 1)
            new_labels[u] = f"h{h:016x}"
        labels = new_labels

    final = sorted(labels[u] for u in nodes)
    s2 = ";".join(final)
    h2 = 2166136261
    for ch in s2.encode("utf-8"):
        h2 ^= ch
        h2 = (h2 * 16777619) & 0xFFFFFFFF
    return f"wl{h2:08x}"


def formula(atoms: Atoms, nodes: list[int]) -> str:
    if len(nodes) == 2:
        syms2 = sorted(atoms[i].symbol for i in nodes)
        if syms2 == ["H", "O"]:
            return "OH"

    syms = [atoms[i].symbol for i in nodes]
    c = Counter(syms)
    keys = []
    if "C" in c:
        keys.append("C")
    if "H" in c:
        keys.append("H")
    for k in sorted(c.keys()):
        if k not in ("C", "H"):
            keys.append(k)
    out = []
    for k in keys:
        n = c[k]
        out.append(k if n == 1 else f"{k}{n}")
    return "".join(out)


def normalize_fragment_smiles(
    atoms: Atoms,
    nodes: list[int],
    cov_edges: set[tuple[int, int]],
    smiles: str,
) -> str:
    if not smiles:
        return smiles

    if len(nodes) == 2:
        node_set = set(nodes)
        n_cov = sum(1 for i, j in cov_edges if i in node_set and j in node_set)
        syms2 = sorted(atoms[i].symbol for i in nodes)
        if syms2 == ["H", "O"] and n_cov == 1:
            return "[OH]"

    return smiles


@dataclass(frozen=True)
class SurfaceSignature:
    n_ads_bonds: int
    slab_coord_hist: tuple[int, ...]
    slab_elements: tuple[str, ...]
    ads_pairs: tuple[tuple[int, int], ...]
    ads_pair_labels: tuple[str, ...]


def surface_signature(
    atoms: Atoms,
    comp_nodes: list[int],
    slab_mask: np.ndarray,
    ads_edges: set[tuple[int, int]],
) -> SurfaceSignature | None:
    comp = set(comp_nodes)
    syms = atoms.get_chemical_symbols()
    per_atom: dict[int, int] = defaultdict(int)
    touched: list[str] = []
    ads_pairs: list[tuple[int, int]] = []
    n_ads = 0

    for i, j in ads_edges:
        if i in comp and slab_mask[j]:
            n_ads += 1
            per_atom[i] += 1
            touched.append(syms[j])
            ads_pairs.append((int(i), int(j)))
        elif j in comp and slab_mask[i]:
            n_ads += 1
            per_atom[j] += 1
            touched.append(syms[i])
            ads_pairs.append((int(j), int(i)))

    if n_ads == 0:
        return None

    sorted_ads_pairs = tuple(sorted(ads_pairs, key=lambda p: (p[0], p[1])))
    ads_pair_labels = tuple(f"{syms[ads]}{ads}-{syms[slab]}{slab}" for ads, slab in sorted_ads_pairs)

    return SurfaceSignature(
        n_ads_bonds=n_ads,
        slab_coord_hist=tuple(sorted(per_atom.values())),
        slab_elements=tuple(sorted(touched)),
        ads_pairs=sorted_ads_pairs,
        ads_pair_labels=ads_pair_labels,
    )


def _single_atom_smiles(atoms: Atoms, nodes: list[int]) -> str:
    if len(nodes) != 1:
        return ""
    sym = atoms[nodes[0]].symbol
    return f"[{sym}]"


def _fallback_atom_list_smiles(atoms: Atoms, nodes: list[int]) -> str:
    zs = [(int(atoms.numbers[a]), a) for a in nodes]
    zs.sort(key=lambda x: (x[0], x[1]))
    syms = atoms.get_chemical_symbols()
    return "".join(f"[{syms[a]}]" for _z, a in zs)


def component_smiles_rdkit_3d(
    atoms: Atoms,
    nodes: list[int],
    *,
    hide_hs: bool = False,
) -> str:
    """
    RDKit geometry-based SMILES (3D bond perception).
    Returns "" on failure or uninformative SMILES.
    """
    smi1 = _single_atom_smiles(atoms, nodes)
    if smi1:
        return smi1

    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return ""

    try:
        rw = Chem.RWMol()
        for a in nodes:
            rw.AddAtom(Chem.Atom(int(atoms.numbers[a])))

        m = rw.GetMol()
        conf = Chem.Conformer(len(nodes))
        for i, a in enumerate(nodes):
            x, y, z = atoms.positions[a]
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        m.AddConformer(conf, assignId=True)

        rdDetermineBonds.DetermineBonds(m)
        Chem.SanitizeMol(m)

        mol_use = m
        if hide_hs and mol_use.GetNumHeavyAtoms() > 0:
            try:
                m2 = Chem.RemoveHs(mol_use)
                if m2.GetNumAtoms() > 0:
                    mol_use = m2
            except Exception:
                pass

        smi = Chem.MolToSmiles(mol_use, canonical=True, allHsExplicit=not hide_hs)

        if smi and smi not in ("*", "[*]"):
            return smi
        return ""
    except Exception:
        return ""


def component_smiles_openbabel_3d(
    atoms: Atoms,
    nodes: list[int],
) -> str:
    """
    OpenBabel geometry-based SMILES (3D bond perception).
    Returns "" on failure or uninformative SMILES.
    """
    smi1 = _single_atom_smiles(atoms, nodes)
    if smi1:
        return smi1

    try:
        from openbabel import openbabel as ob
    except Exception:
        return ""

    try:
        mol = ob.OBMol()
        mol.BeginModify()

        for a in nodes:
            at = mol.NewAtom()
            at.SetAtomicNum(int(atoms.numbers[a]))
            x, y, z = atoms.positions[a]
            at.SetVector(float(x), float(y), float(z))

        mol.EndModify()

        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        ob.OBAromTyper().AssignAromaticFlags(mol)

        conv = ob.OBConversion()
        conv.SetOutFormat("can")
        smi = conv.WriteString(mol).strip()
        smi = smi.split()[0].strip() if smi else ""

        if smi and smi not in ("*", "[*]"):
            return smi
        return ""
    except Exception:
        return ""


def component_smiles_best_effort(
    atoms: Atoms,
    nodes: list[int],
) -> str:
    """
    Best-effort, chemistry-first SMILES from 3D geometry.
    - Prefer RDKit 3D
    - Fallback: OpenBabel 3D
    - Final fallback: concatenated atomic SMILES like [H][H][O]
    """
    smi = component_smiles_rdkit_3d(atoms, nodes)
    if smi:
        return smi

    smi = component_smiles_openbabel_3d(atoms, nodes)
    if smi:
        return smi

    return _fallback_atom_list_smiles(atoms, nodes)


def component_smiles_from_edges_rdkit(
    atoms: Atoms,
    nodes: list[int],
    cov_edges: set[tuple[int, int]],
    *,
    sanitize: bool = True,
    canonical: bool = True,
    hide_hs: bool = False,
) -> tuple[str | None, bool]:
    """
    Build RDKit Mol using confirmed covalent edges (topology-first),
    then output SMILES. All covalent edges are treated as SINGLE bonds.

    Returns: (smiles, failed)
      - smiles: canonical SMILES string or None
      - failed: True if RDKit was available but conversion/sanitization failed
    """
    smi1 = _single_atom_smiles(atoms, nodes)
    if smi1:
        return smi1, False

    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return None, False

    try:
        node_set = set(nodes)
        idx_map: dict[int, int] = {}
        rw = Chem.RWMol()
        for a in nodes:
            rid = rw.AddAtom(Chem.Atom(int(atoms.numbers[a])))
            idx_map[a] = rid

        for i, j in cov_edges:
            if i in node_set and j in node_set:
                ri = idx_map[i]
                rj = idx_map[j]
                if rw.GetBondBetweenAtoms(ri, rj) is None:
                    rw.AddBond(ri, rj, Chem.BondType.SINGLE)

        mol = rw.GetMol()
        if sanitize:
            Chem.SanitizeMol(mol)

        mol_use = mol
        if hide_hs and mol_use.GetNumHeavyAtoms() > 0:
            try:
                m2 = Chem.RemoveHs(mol_use)
                if m2.GetNumAtoms() > 0:
                    mol_use = m2
            except Exception:
                pass

        smi = Chem.MolToSmiles(mol_use, canonical=canonical, allHsExplicit=not hide_hs)
        if smi and smi not in ("*", "[*]"):
            return smi, False
        return None, True
    except Exception:
        return None, True


def is_suspicious_smiles(
    smi: str,
    *,
    allow_charged: bool = False,
    allow_dot: bool = False,
) -> bool:
    if not smi:
        return True
    if (not allow_dot) and ("." in smi):
        return True
    if (not allow_charged) and (("+" in smi) or ("-" in smi)):
        return True
    return False


class SmilesStrategy(ABC):
    @abstractmethod
    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        ...


@dataclass
class SmilesRDKit3DStrategy(SmilesStrategy):
    hide_hs: bool = False

    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        return component_smiles_rdkit_3d(atoms, nodes, hide_hs=self.hide_hs)


@dataclass
class SmilesOpenBabel3DStrategy(SmilesStrategy):
    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        return component_smiles_openbabel_3d(atoms, nodes)


@dataclass
class SmilesRDKitTopologyStrategy(SmilesStrategy):
    sanitize: bool = True
    hide_hs: bool = False

    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        smi, _failed = component_smiles_from_edges_rdkit(
            atoms, nodes, cov_edges, sanitize=self.sanitize, canonical=True, hide_hs=self.hide_hs
        )
        return smi or ""


@dataclass
class SmilesComboStrategy(SmilesStrategy):
    strategies: list[SmilesStrategy]
    treat_suspicious_as_failure: bool = True
    allow_charged: bool | None = None
    allow_dot: bool | None = None
    prefer_charged: bool = False

    def _effective_flags(self) -> tuple[bool, bool]:
        if self.allow_charged is None:
            allow_charged = not self.treat_suspicious_as_failure
        else:
            allow_charged = self.allow_charged
        if self.allow_dot is None:
            allow_dot = not self.treat_suspicious_as_failure
        else:
            allow_dot = self.allow_dot
        return allow_charged, allow_dot

    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        allow_charged, allow_dot = self._effective_flags()

        first_ok = ""
        for strat in self.strategies:
            smi = strat.compute(atoms, nodes, cov_edges) or ""
            if not smi:
                continue

            if is_suspicious_smiles(smi, allow_charged=allow_charged, allow_dot=allow_dot):
                continue

            is_charged = ("+" in smi) or ("-" in smi)

            if self.prefer_charged and is_charged and allow_charged:
                return smi

            if not first_ok:
                first_ok = smi

            if not self.prefer_charged:
                return smi

        return first_ok


@dataclass
class SmilesBestEffortStrategy(SmilesStrategy):
    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        return component_smiles_best_effort(atoms, nodes)


@dataclass
class SmilesEdgesRDKitStrategy(SmilesStrategy):
    sanitize: bool = True
    hide_hs: bool = False
    fallback: SmilesStrategy | None = None

    def compute(self, atoms: Atoms, nodes: list[int], cov_edges: set[tuple[int, int]]) -> str:
        smi = SmilesRDKitTopologyStrategy(sanitize=self.sanitize, hide_hs=self.hide_hs).compute(
            atoms, nodes, cov_edges
        )
        if smi:
            return smi
        if self.fallback is not None:
            return self.fallback.compute(atoms, nodes, cov_edges)
        return ""
