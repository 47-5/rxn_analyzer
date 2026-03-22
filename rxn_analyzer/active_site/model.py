from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field


def _normalize_index(raw: object, *, field_name: str, index_base: int) -> int:
    try:
        v = int(raw)
    except Exception as e:
        raise ValueError(f"{field_name} has non-integer atom index: {raw!r}") from e

    if index_base == 0:
        if v < 0:
            raise ValueError(f"{field_name} has invalid 0-based index {v} (must be >= 0).")
        return v

    if index_base == 1:
        if v < 1:
            raise ValueError(f"{field_name} has invalid 1-based index {v} (must be >= 1).")
        return v - 1

    raise ValueError(f"index_base must be 0 or 1, got {index_base}.")


def _normalize_indices(raw: object, *, field_name: str, index_base: int) -> frozenset[int]:
    if raw is None:
        return frozenset()
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Iterable):
        raise ValueError(f"{field_name} must be an iterable of integers.")

    out: list[int] = []
    for item in raw:
        out.append(_normalize_index(item, field_name=field_name, index_base=index_base))
    return frozenset(out)


def _normalize_str_set(raw: object, *, field_name: str) -> frozenset[str]:
    if raw is None:
        return frozenset()
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Iterable):
        raise ValueError(f"{field_name} must be an iterable of strings.")
    return frozenset(str(x).strip() for x in raw if str(x).strip())


@dataclass(frozen=True)
class ActiveSite:
    site_id: str
    site_family: str
    core_members: frozenset[int]
    initial_state_members: frozenset[int] = frozenset()
    allowed_state_elements: frozenset[str] = frozenset()
    max_state_members: int | None = None
    rule_profile: str = "default"
    metadata: dict[str, object] = field(default_factory=dict)

    @staticmethod
    def from_dict(d: Mapping[str, object], *, index_base: int = 0) -> "ActiveSite":
        site_id = str(d.get("id", d.get("site_id", ""))).strip()
        if not site_id:
            raise ValueError(f"Reactive site missing id/site_id: {d}")

        family = str(d.get("family", d.get("type", d.get("site_family", "")))).strip()
        if not family:
            raise ValueError(f"Reactive site '{site_id}' missing family/type.")

        core = _normalize_indices(d.get("core_members", d.get("core_atoms", [])),
                                  field_name=f"ActiveSite '{site_id}' core_members",
                                  index_base=index_base)
        if not core:
            raise ValueError(f"Reactive site '{site_id}' has empty core_members.")

        initial = _normalize_indices(
            d.get("initial_state_members", d.get("state_members", [])),
            field_name=f"ActiveSite '{site_id}' initial_state_members",
            index_base=index_base,
        )
        allowed = _normalize_str_set(
            d.get("allowed_state_elements", []),
            field_name=f"ActiveSite '{site_id}' allowed_state_elements",
        )

        max_state_members_raw = d.get("max_state_members", None)
        max_state_members = None if max_state_members_raw is None else int(max_state_members_raw)
        if max_state_members is not None and max_state_members < 0:
            raise ValueError(f"Active site '{site_id}' max_state_members must be >= 0.")

        rule_profile = str(d.get("rule_profile", "default")).strip() or "default"
        metadata = dict(d.get("metadata", {})) if isinstance(d.get("metadata", {}), Mapping) else {}

        return ActiveSite(
            site_id=site_id,
            site_family=family,
            core_members=core,
            initial_state_members=initial,
            allowed_state_elements=allowed,
            max_state_members=max_state_members,
            rule_profile=rule_profile,
            metadata=metadata,
        )


@dataclass(frozen=True)
class ActiveSiteDefinition:
    sites: tuple[ActiveSite, ...] = ()
    strict_core_validation: bool = True

    def __post_init__(self) -> None:
        ids = [s.site_id for s in self.sites]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate reactive site ids found: {ids}")

    @classmethod
    def from_mapping(
        cls,
        cfg: Mapping[str, object],
        *,
        index_base: int = 0,
        strict_core_validation: bool = True,
    ) -> "ActiveSiteDefinition":
        raw = cfg.get("sites", cfg.get("active_sites", cfg.get("reactive_sites", [])))
        if raw is None:
            raw = []
        if not isinstance(raw, list):
            raise ValueError("Reactive site config must contain a list under key 'sites'.")
        sites = tuple(ActiveSite.from_dict(item, index_base=index_base) for item in raw)
        return cls(sites=sites, strict_core_validation=bool(strict_core_validation))

    @classmethod
    def from_yaml(
        cls,
        path: str,
        *,
        index_base: int = 0,
        strict_core_validation: bool = True,
    ) -> "ActiveSiteDefinition":
        try:
            import yaml  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "Loading YAML reactive-site config requires PyYAML. Please `pip install pyyaml`."
            ) from e

        with open(path, "r", encoding="utf-8") as f:
            obj = yaml.safe_load(f) or {}
        if not isinstance(obj, Mapping):
            raise ValueError("Reactive site YAML root must be a mapping/object.")
        return cls.from_mapping(
            obj,
            index_base=index_base,
            strict_core_validation=strict_core_validation,
        )


@dataclass(frozen=True)
class ActiveSiteStateFrame:
    frame: int
    site_id: str
    site_family: str
    core_members: tuple[int, ...]
    intrinsic_members: tuple[int, ...]
    incorporated_members: tuple[int, ...]
    associated_members: tuple[int, ...]
    state_members: tuple[int, ...]
    intrinsic_formula: str
    intrinsic_topology: str
    incorporated_formula: str
    incorporated_topology: str
    associated_formula: str
    associated_topology: str
    state_formula: str
    state_topology: str
    state_label: str
    associated_species_labels: tuple[str, ...]
    attached_components: tuple[tuple[int, ...], ...]
    descriptors: dict[str, object]


@dataclass(frozen=True)
class ActiveSiteEvent:
    event_id: int
    frame: int
    site_id: str
    event_type: str
    old_state: str
    new_state: str
    added_members: tuple[int, ...]
    removed_members: tuple[int, ...]
    associated_species_before: tuple[str, ...]
    associated_species_after: tuple[str, ...]
    evidence_bonds: list[dict]


@dataclass(frozen=True)
class JointActiveSiteReaction:
    event_id: int
    frame: int
    site_id: str
    site_state_before: str
    site_state_after: str
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    associated_species_before: tuple[str, ...]
    associated_species_after: tuple[str, ...]
    evidence_bonds: list[dict]


@dataclass(frozen=True)
class SiteReactionCouplingRow:
    frame: int
    site_id: str
    site_family: str
    site_state_before: str
    site_state_after: str
    site_changed: bool
    has_transform: bool
    transform_event_ids: tuple[int, ...]
    transform_types: tuple[str, ...]
    transform_reactants: tuple[str, ...]
    transform_products: tuple[str, ...]
    associated_species_before: tuple[str, ...]
    associated_species_after: tuple[str, ...]
    coupling_type: str
    link_strength: str
    evidence_bonds: list[dict]


# Backward aliases for older internal names.
ReactiveSite = ActiveSite
ReactiveSiteDefinition = ActiveSiteDefinition
ReactiveSiteStateFrame = ActiveSiteStateFrame
ReactiveSiteEvent = ActiveSiteEvent
JointReactiveSiteReaction = JointActiveSiteReaction
