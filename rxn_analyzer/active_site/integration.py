from __future__ import annotations

from .model import ActiveSiteStateFrame


def build_active_site_memberships(
    comps: list[list[int]],
    states: list[ActiveSiteStateFrame],
) -> list[list[dict]]:
    memberships: list[list[dict]] = []
    for comp in comps:
        comp_set = set(comp)
        comp_items: list[dict] = []
        for state in states:
            overlap_rows: list[dict] = []
            for role_name, atoms in (
                ("intrinsic", set(state.intrinsic_members)),
                ("incorporated", set(state.incorporated_members)),
                ("associated", set(state.associated_members)),
            ):
                overlap = sorted(comp_set & atoms)
                if overlap:
                    overlap_rows.append({"role": role_name, "overlap_atoms": overlap})
            if overlap_rows:
                comp_items.append(
                    {
                        "site_id": state.site_id,
                        "site_family": state.site_family,
                        "memberships": overlap_rows,
                    }
                )
        memberships.append(comp_items)
    return memberships


def summarize_active_site_memberships(
    memberships: list[list[dict]],
) -> tuple[list[str], list[list[str]], list[list[str]], list[bool], list[bool]]:
    role_priority = {"intrinsic": 0, "incorporated": 1, "associated": 2}
    summary_roles: list[str] = []
    summary_site_ids: list[list[str]] = []
    summary_labels: list[list[str]] = []
    is_site_owned: list[bool] = []
    is_site_associated: list[bool] = []

    for comp_items in memberships:
        roles_found: set[str] = set()
        site_ids: list[str] = []
        labels: list[str] = []

        for item in comp_items:
            site_id = str(item.get("site_id", ""))
            site_ids.append(site_id)
            memberships_list = item.get("memberships", [])
            roles = sorted(
                {
                    str(entry.get("role", ""))
                    for entry in memberships_list
                    if str(entry.get("role", ""))
                },
                key=lambda x: (role_priority.get(x, 99), x),
            )
            roles_found.update(roles)
            if site_id and roles:
                labels.append(f"{site_id}:{'+'.join(roles)}")
            elif site_id:
                labels.append(site_id)

        normalized_site_ids = sorted({x for x in site_ids if x})
        normalized_labels = sorted(set(labels))

        owned = any(r in {"intrinsic", "incorporated"} for r in roles_found)
        associated = "associated" in roles_found

        if not roles_found:
            role_text = "none"
        elif len(roles_found) == 1:
            role_text = next(iter(roles_found))
        else:
            role_text = "mixed"

        summary_roles.append(role_text)
        summary_site_ids.append(normalized_site_ids)
        summary_labels.append(normalized_labels)
        is_site_owned.append(owned)
        is_site_associated.append(associated)

    return summary_roles, summary_site_ids, summary_labels, is_site_owned, is_site_associated
