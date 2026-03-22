from __future__ import annotations

from dataclasses import dataclass

from .model import ActiveSiteStateFrame


@dataclass(frozen=True)
class ActiveSiteFrameBatch:
    states: list[ActiveSiteStateFrame]


@dataclass(frozen=True)
class MembershipCandidates:
    valid_core: frozenset[int]
    direct_contacts: set[int]
    intrinsic_seed_members: set[int]
    intrinsic_members: set[int]
    candidate_component_atoms: dict[tuple[int, ...], set[int]]


@dataclass(frozen=True)
class DynamicMembership:
    incorporated_members: set[int]
    associated_members: set[int]
    remaining_component_keys: list[tuple[int, ...]]


@dataclass(frozen=True)
class LimitedMembership:
    intrinsic_members: set[int]
    incorporated_members: set[int]
    associated_members: set[int]
    state_members: set[int]


@dataclass(frozen=True)
class AssociatedComponentInfo:
    component_keys: list[tuple[int, ...]]
    species_labels: list[str]


@dataclass(frozen=True)
class ActiveSiteFrameContext:
    memberships: list[list[dict]]
    roles: list[str]
    site_ids: list[list[str]]
    summaries: list[list[str]]
    is_site_owned: list[bool]
    is_site_associated: list[bool]
