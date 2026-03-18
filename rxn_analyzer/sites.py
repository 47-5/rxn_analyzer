from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from collections.abc import Iterable, Mapping
import numpy as np


def _normalize_index(raw: object, *, site_id: str, index_base: int) -> int:
    """
    Normalize site atom index to 0-based internal index.
    """
    try:
        v = int(raw)
    except Exception as e:
        raise ValueError(f"Site '{site_id}' has non-integer atom index: {raw!r}") from e

    if index_base == 0:
        if v < 0:
            raise ValueError(
                f"Site '{site_id}' has invalid 0-based index {v} (must be >= 0)."
            )
        return v

    if index_base == 1:
        if v < 1:
            raise ValueError(
                f"Site '{site_id}' has invalid 1-based index {v} (must be >= 1)."
            )
        return v - 1

    raise ValueError(f"index_base must be 0 or 1, got {index_base}.")


@dataclass(frozen=True)
class Site:
    site_id: str
    site_type: str
    atom_indices: frozenset[int]
    name: str = ""
    priority: int = 0

    @staticmethod
    def from_dict(d: Mapping[str, object], *, index_base: int = 0) -> "Site":
        sid = str(d.get("id", d.get("site_id", ""))).strip()
        stype = str(d.get("type", d.get("site_type", ""))).strip()
        if not sid:
            raise ValueError(f"Site missing id/site_id: {d}")
        if not stype:
            raise ValueError(f"Site missing type/site_type: {d}")

        atoms_raw = d.get("atoms", d.get("indices", []))
        if atoms_raw is None:
            atoms_raw = []

        if isinstance(atoms_raw, (str, bytes)) or not isinstance(atoms_raw, Iterable):
            raise ValueError(f"Site '{sid}' atoms/indices must be an iterable of integers.")

        converted: list[int] = []
        for x in atoms_raw:
            converted.append(_normalize_index(x, site_id=sid, index_base=index_base))

        if len(converted) == 0:
            raise ValueError(f"Site '{sid}' has empty atoms/indices.")

        seen = set()
        dups = set()
        for a in converted:
            if a in seen:
                dups.add(a)
            seen.add(a)
        if dups:
            raise ValueError(
                f"Site '{sid}' has duplicate atom indices after base conversion: {sorted(dups)}"
            )

        name = str(d.get("name", "") or "")
        priority = int(d.get("priority", 0))

        return Site(
            site_id=sid,
            site_type=stype,
            atom_indices=frozenset(converted),
            name=name,
            priority=priority,
        )


@dataclass(frozen=True)
class SiteAssignment:
    primary_site_id: str
    primary_site_type: str
    # only top-score ties, sorted by site_id for determinism
    candidate_sites: tuple[tuple[str, str], ...]
    ambiguous: bool
    n_ads_bonds: int
    touched_host_atoms: tuple[int, ...]


@dataclass(frozen=True)
class SiteDefinition:
    sites: tuple[Site, ...] = ()
    strict_index_validation: bool = True

    def __post_init__(self):
        ids = [s.site_id for s in self.sites]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate site_id found: {ids}")

    @classmethod
    def from_mapping(
        cls,
        cfg: Mapping[str, object],
        *,
        index_base: int = 0,                # 0 or 1
        strict_index_validation: bool = True,
    ) -> "SiteDefinition":
        if index_base not in (0, 1):
            raise ValueError(f"index_base must be 0 or 1, got {index_base}.")

        raw = cfg.get("sites", cfg.get("site_definitions", []))
        if raw is None:
            raw = []
        if not isinstance(raw, list):
            raise ValueError("Site config must contain a list under key 'sites'.")

        sites: list[Site] = []
        for item in raw:
            if not isinstance(item, Mapping):
                raise ValueError(f"Each site entry must be a mapping, got: {item!r}")
            sites.append(Site.from_dict(item, index_base=index_base))

        return cls(
            sites=tuple(sites),
            strict_index_validation=bool(strict_index_validation),
        )

    @classmethod
    def from_yaml(
        cls,
        path: str,
        *,
        index_base: int = 0,                # 0 or 1
        strict_index_validation: bool = True,
    ) -> "SiteDefinition":
        try:
            import yaml  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "Loading YAML site config requires PyYAML. Please `pip install pyyaml`."
            ) from e

        with open(path, "r", encoding="utf-8") as f:
            obj = yaml.safe_load(f) or {}
        if not isinstance(obj, Mapping):
            raise ValueError("YAML root must be a mapping/object.")
        return cls.from_mapping(
            obj,
            index_base=index_base,
            strict_index_validation=strict_index_validation,
        )

    def _build_atom_to_sites(self, slab_mask: np.ndarray) -> dict[int, list[Site]]:
        """
        Build host-atom -> site list map.
        Strict mode:
          - out-of-range indices -> raise
          - indices not marked as host/slab -> raise
        Non-strict mode:
          - invalid ones are skipped
        """
        n = int(len(slab_mask))
        atom_to_sites: dict[int, list[Site]] = defaultdict(list)

        idx_errors: dict[str, list[int]] = defaultdict(list)
        nonhost_errors: dict[str, list[int]] = defaultdict(list)

        for s in self.sites:
            for a in s.atom_indices:
                if a < 0 or a >= n:
                    idx_errors[s.site_id].append(a)
                    continue
                if not bool(slab_mask[a]):
                    nonhost_errors[s.site_id].append(a)
                    continue
                atom_to_sites[a].append(s)

        if (idx_errors or nonhost_errors) and self.strict_index_validation:
            lines: list[str] = ["Invalid site definition for current frame:"]
            if idx_errors:
                lines.append("  Out-of-range indices:")
                for sid in sorted(idx_errors):
                    bad = sorted(set(idx_errors[sid]))
                    lines.append(f"    - {sid}: {bad} (valid range: 0..{max(0, n - 1)})")
            if nonhost_errors:
                lines.append("  Indices not in host/slab mask:")
                for sid in sorted(nonhost_errors):
                    bad = sorted(set(nonhost_errors[sid]))
                    lines.append(f"    - {sid}: {bad}")
            raise ValueError("\n".join(lines))

        return atom_to_sites

    def assign_component(
        self,
        comp_nodes: list[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
    ) -> SiteAssignment | None:
        if not self.sites:
            return None

        comp = set(comp_nodes)
        if not comp:
            return None

        atom_to_sites = self._build_atom_to_sites(slab_mask)
        if not atom_to_sites:
            return None

        by_id = {s.site_id: s for s in self.sites}
        n = int(len(slab_mask))

        bond_hits: dict[str, int] = defaultdict(int)
        site_host_atoms: dict[str, set[int]] = defaultdict(set)
        touched: list[int] = []
        n_ads = 0

        for i, j in ads_edges:
            host = None
            if i in comp and 0 <= j < n and bool(slab_mask[j]):
                host = j
            elif j in comp and 0 <= i < n and bool(slab_mask[i]):
                host = i

            if host is None:
                continue

            n_ads += 1
            touched.append(host)

            for s in atom_to_sites.get(host, []):
                bond_hits[s.site_id] += 1
                site_host_atoms[s.site_id].add(host)

        if n_ads == 0 or not bond_hits:
            return None

        # score: n_ads_bond_hits > n_unique_host_atoms > site.priority
        scored: list[tuple[str, int, int, int]] = []
        for sid, bh in bond_hits.items():
            ah = len(site_host_atoms.get(sid, set()))
            pr = int(by_id[sid].priority)
            scored.append((sid, int(bh), int(ah), pr))

        best_key = max((bh, ah, pr) for _sid, bh, ah, pr in scored)
        winners = [(sid, bh, ah, pr) for sid, bh, ah, pr in scored if (bh, ah, pr) == best_key]
        winners.sort(key=lambda x: x[0])

        primary_sid = winners[0][0]
        primary_stype = by_id[primary_sid].site_type
        candidate_sites = tuple((sid, by_id[sid].site_type) for sid, *_ in winners)

        return SiteAssignment(
            primary_site_id=primary_sid,
            primary_site_type=primary_stype,
            candidate_sites=candidate_sites,
            ambiguous=(len(winners) > 1),
            n_ads_bonds=n_ads,
            touched_host_atoms=tuple(sorted(set(touched))),
        )