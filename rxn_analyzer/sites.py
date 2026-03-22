from __future__ import annotations
from dataclasses import dataclass, field
from itertools import combinations
from collections import defaultdict
from collections.abc import Iterable, Mapping
import numpy as np
from ase import Atoms


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


# Preferred semantic aliases for static geometric adsorption sites.
GeometricSite = Site


def _normalize_str_tuple(raw: object, *, field_name: str) -> tuple[str, ...]:
    if raw is None:
        return ()
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Iterable):
        raise ValueError(f"{field_name} must be an iterable of strings.")

    out: list[str] = []
    for item in raw:
        token = str(item).strip().lower()
        if token:
            out.append(token)
    return tuple(out)


@dataclass(frozen=True)
class AutoSiteConfig:
    enabled: bool = False
    allowed_types: tuple[str, ...] = ("top", "bridge", "fcc", "hcp", "fourfold")
    bridge_requires_edge: bool = True
    hollow_requires_complete_triangle: bool = True
    fourfold_requires_cycle: bool = True
    layer_tolerance: float = 0.75
    subsurface_match_tolerance: float = 1.25
    type_priorities: dict[str, int] = field(
        default_factory=lambda: {
            "top": 10,
            "bridge": 20,
            "hollow": 30,
            "fcc": 30,
            "hcp": 30,
            "fourfold": 40,
        }
    )

    @classmethod
    def from_mapping(cls, cfg: Mapping[str, object] | None) -> "AutoSiteConfig":
        if cfg is None:
            return cls()

        enabled = bool(cfg.get("enabled", True))
        allowed_types = _normalize_str_tuple(
            cfg.get("allowed_types", ("top", "bridge", "fcc", "hcp", "fourfold")),
            field_name="auto_site.allowed_types",
        )

        priorities_raw = cfg.get("type_priorities", {})
        priorities = {
            "top": 10,
            "bridge": 20,
            "hollow": 30,
            "fcc": 30,
            "hcp": 30,
            "fourfold": 40,
        }
        if isinstance(priorities_raw, Mapping):
            for key, value in priorities_raw.items():
                priorities[str(key).strip().lower()] = int(value)

        return cls(
            enabled=enabled,
            allowed_types=allowed_types or ("top", "bridge", "fcc", "hcp", "fourfold"),
            bridge_requires_edge=bool(cfg.get("bridge_requires_edge", True)),
            hollow_requires_complete_triangle=bool(cfg.get("hollow_requires_complete_triangle", True)),
            fourfold_requires_cycle=bool(cfg.get("fourfold_requires_cycle", True)),
            layer_tolerance=float(cfg.get("layer_tolerance", 0.75)),
            subsurface_match_tolerance=float(cfg.get("subsurface_match_tolerance", 1.25)),
            type_priorities=priorities,
        )

    def priority_for(self, site_type: str) -> int:
        return int(self.type_priorities.get(site_type.lower().strip(), 0))


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
    auto_site_config: AutoSiteConfig = field(default_factory=AutoSiteConfig)

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

        raw = cfg.get("sites", cfg.get("geometric_sites", cfg.get("site_definitions", [])))
        if raw is None:
            raw = []
        if not isinstance(raw, list):
            raise ValueError("Site config must contain a list under key 'sites'.")

        sites: list[Site] = []
        for item in raw:
            if not isinstance(item, Mapping):
                raise ValueError(f"Each site entry must be a mapping, got: {item!r}")
            sites.append(Site.from_dict(item, index_base=index_base))

        auto_cfg_raw = cfg.get("auto_site", cfg.get("auto_detection", cfg.get("auto_generate", None)))
        if isinstance(auto_cfg_raw, bool):
            auto_cfg = AutoSiteConfig(enabled=bool(auto_cfg_raw))
        elif isinstance(auto_cfg_raw, Mapping):
            auto_cfg = AutoSiteConfig.from_mapping(auto_cfg_raw)
        elif auto_cfg_raw is None:
            auto_cfg = AutoSiteConfig()
        else:
            raise ValueError("auto_site/auto_detection/auto_generate must be a mapping or boolean.")

        return cls(
            sites=tuple(sites),
            strict_index_validation=bool(strict_index_validation),
            auto_site_config=auto_cfg,
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

    @staticmethod
    def _sorted_edge(i: int, j: int) -> tuple[int, int]:
        return (i, j) if i <= j else (j, i)

    @staticmethod
    def _ads_host_contacts(
        comp: set[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
    ) -> tuple[dict[int, int], int]:
        counts: dict[int, int] = defaultdict(int)
        total = 0
        n = int(len(slab_mask))

        for i, j in ads_edges:
            host = None
            if i in comp and 0 <= j < n and bool(slab_mask[j]):
                host = j
            elif j in comp and 0 <= i < n and bool(slab_mask[i]):
                host = i

            if host is None:
                continue

            total += 1
            counts[host] += 1

        return counts, total

    @staticmethod
    def _surface_normal(atoms: Atoms, slab_mask: np.ndarray, comp: set[int]) -> np.ndarray | None:
        slab_indices = np.flatnonzero(slab_mask)
        if len(slab_indices) < 3:
            return None

        slab_pos = atoms.positions[slab_indices]
        centered = slab_pos - slab_pos.mean(axis=0)
        cov = centered.T @ centered
        vals, vecs = np.linalg.eigh(cov)
        normal = np.array(vecs[:, int(np.argmin(vals))], dtype=float)
        norm = float(np.linalg.norm(normal))
        if norm <= 0.0:
            return None
        normal /= norm

        if comp:
            ads_center = atoms.positions[sorted(comp)].mean(axis=0)
            slab_center = slab_pos.mean(axis=0)
            if float(np.dot(ads_center - slab_center, normal)) < 0.0:
                normal *= -1.0

        return normal

    @staticmethod
    def _layer_groups(
        atoms: Atoms,
        slab_mask: np.ndarray,
        normal: np.ndarray,
        layer_tolerance: float,
    ) -> tuple[list[list[int]], dict[int, float]]:
        slab_indices = [int(i) for i in np.flatnonzero(slab_mask)]
        heights = {i: float(np.dot(atoms.positions[i], normal)) for i in slab_indices}
        ordered = sorted(slab_indices, key=lambda idx: heights[idx], reverse=True)

        layers: list[list[int]] = []
        for idx in ordered:
            if not layers:
                layers.append([idx])
                continue

            ref_h = float(np.mean([heights[x] for x in layers[-1]]))
            if abs(heights[idx] - ref_h) <= layer_tolerance:
                layers[-1].append(idx)
            else:
                layers.append([idx])

        return layers, heights

    @staticmethod
    def _site_id(site_type: str, atom_indices: tuple[int, ...]) -> str:
        joined = "_".join(str(i) for i in atom_indices)
        return f"auto_{site_type}_{joined}"

    def _classify_hollow_type(
        self,
        atoms: Atoms,
        slab_mask: np.ndarray,
        comp: set[int],
        atom_indices: tuple[int, int, int],
    ) -> str:
        allowed = set(self.auto_site_config.allowed_types)
        if "fcc" not in allowed and "hcp" not in allowed:
            return "hollow"

        normal = self._surface_normal(atoms, slab_mask, comp)
        if normal is None:
            return "hollow"

        layers, heights = self._layer_groups(
            atoms,
            slab_mask,
            normal,
            self.auto_site_config.layer_tolerance,
        )
        if not layers:
            return "hollow"

        atom_set = set(atom_indices)
        host_layer_idx = None
        for idx, layer in enumerate(layers):
            if atom_set.issubset(layer):
                host_layer_idx = idx
                break

        if host_layer_idx is None:
            triplet_heights = [heights[a] for a in atom_indices]
            if max(triplet_heights) - min(triplet_heights) > self.auto_site_config.layer_tolerance:
                return "hollow"
            host_mean = float(np.mean(triplet_heights))
            host_layer_idx = min(
                range(len(layers)),
                key=lambda idx: abs(float(np.mean([heights[a] for a in layers[idx]])) - host_mean),
            )

        if host_layer_idx + 1 >= len(layers):
            return "fcc" if "fcc" in allowed else "hollow"


        center = atoms.positions[list(atom_indices)].mean(axis=0)
        subsurface = layers[host_layer_idx + 1]
        for idx in subsurface:
            delta = atoms.positions[idx] - center
            lateral = delta - normal * float(np.dot(delta, normal))
            if float(np.linalg.norm(lateral)) <= self.auto_site_config.subsurface_match_tolerance:
                return "hcp" if "hcp" in allowed else "hollow"

        return "fcc" if "fcc" in allowed else "hollow"

    def _build_auto_sites(
        self,
        atoms: Atoms | None,
        comp: set[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
        slab_edges: set[tuple[int, int]] | None,
    ) -> tuple[Site, ...]:
        if not self.auto_site_config.enabled:
            return ()

        host_counts, _ = self._ads_host_contacts(comp, slab_mask, ads_edges)
        touched_hosts = tuple(sorted(host_counts))
        if not touched_hosts:
            return ()

        allowed = set(self.auto_site_config.allowed_types)
        slab_edge_set = {self._sorted_edge(i, j) for i, j in (slab_edges or set())}
        out: dict[str, Site] = {}

        if "top" in allowed:
            for host in touched_hosts:
                site = Site(
                    site_id=self._site_id("top", (host,)),
                    site_type="top",
                    atom_indices=frozenset((host,)),
                    name=f"auto top at slab atom {host}",
                    priority=self.auto_site_config.priority_for("top"),
                )
                out[site.site_id] = site

        if "bridge" in allowed and len(touched_hosts) >= 2:
            for pair in combinations(touched_hosts, 2):
                if self.auto_site_config.bridge_requires_edge and slab_edges is not None:
                    if self._sorted_edge(*pair) not in slab_edge_set:
                        continue
                site = Site(
                    site_id=self._site_id("bridge", pair),
                    site_type="bridge",
                    atom_indices=frozenset(pair),
                    name=f"auto bridge between slab atoms {pair[0]} and {pair[1]}",
                    priority=self.auto_site_config.priority_for("bridge"),
                )
                out[site.site_id] = site

        hollow_allowed = {"hollow", "fcc", "hcp"} & allowed
        if hollow_allowed and len(touched_hosts) >= 3:
            for trio in combinations(touched_hosts, 3):
                if self.auto_site_config.hollow_requires_complete_triangle and slab_edges is not None:
                    trio_edges = {
                        self._sorted_edge(trio[0], trio[1]),
                        self._sorted_edge(trio[0], trio[2]),
                        self._sorted_edge(trio[1], trio[2]),
                    }
                    if not trio_edges.issubset(slab_edge_set):
                        continue

                site_type = "hollow"
                if atoms is not None:
                    site_type = self._classify_hollow_type(atoms, slab_mask, comp, trio)
                if site_type not in allowed and "hollow" in allowed:
                    site_type = "hollow"
                if site_type not in allowed:
                    continue

                site = Site(
                    site_id=self._site_id(site_type, trio),
                    site_type=site_type,
                    atom_indices=frozenset(trio),
                    name=f"auto {site_type} site on slab atoms {', '.join(str(x) for x in trio)}",
                    priority=self.auto_site_config.priority_for(site_type),
                )
                out[site.site_id] = site

        if "fourfold" in allowed and len(touched_hosts) >= 4:
            for quartet in combinations(touched_hosts, 4):
                if self.auto_site_config.fourfold_requires_cycle and slab_edges is not None:
                    deg = {atom: 0 for atom in quartet}
                    edge_count = 0
                    for i, j in combinations(quartet, 2):
                        if self._sorted_edge(i, j) in slab_edge_set:
                            deg[i] += 1
                            deg[j] += 1
                            edge_count += 1
                    if edge_count < 4 or min(deg.values()) < 2:
                        continue

                site = Site(
                    site_id=self._site_id("fourfold", quartet),
                    site_type="fourfold",
                    atom_indices=frozenset(quartet),
                    name=f"auto fourfold site on slab atoms {', '.join(str(x) for x in quartet)}",
                    priority=self.auto_site_config.priority_for("fourfold"),
                )
                out[site.site_id] = site

        return tuple(sorted(out.values(), key=lambda s: s.site_id))

    def assign_component(
        self,
        comp_nodes: list[int],
        slab_mask: np.ndarray,
        ads_edges: set[tuple[int, int]],
        atoms: Atoms | None = None,
        slab_edges: set[tuple[int, int]] | None = None,
    ) -> SiteAssignment | None:
        if not self.sites and not self.auto_site_config.enabled:
            return None

        comp = set(comp_nodes)
        if not comp:
            return None

        atom_to_sites = self._build_atom_to_sites(slab_mask)
        auto_sites = self._build_auto_sites(atoms, comp, slab_mask, ads_edges, slab_edges)
        for site in auto_sites:
            for atom_idx in site.atom_indices:
                atom_to_sites.setdefault(atom_idx, []).append(site)

        if not atom_to_sites:
            return None

        by_id = {s.site_id: s for s in self.sites}
        by_id.update({s.site_id: s for s in auto_sites})

        host_counts, n_ads = self._ads_host_contacts(comp, slab_mask, ads_edges)
        if n_ads == 0:
            return None

        bond_hits: dict[str, int] = defaultdict(int)
        site_host_atoms: dict[str, set[int]] = defaultdict(set)

        for host, count in host_counts.items():
            for s in atom_to_sites.get(host, []):
                bond_hits[s.site_id] += int(count)
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
            touched_host_atoms=tuple(sorted(host_counts)),
        )


# Preferred semantic aliases for static geometric adsorption sites.
GeometricSiteDefinition = SiteDefinition
GeometricSiteAssignment = SiteAssignment
