"""Configuration loading for the runtime analyzer.

This module only accepts the new config vocabulary:
- `host_definition`
- `geometric_site`
- `active_site`

It is responsible for:
- reading YAML
- applying CLI/runtime overrides
- rejecting legacy keys
- building runtime objects used by `runner.py`
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, fields
from pathlib import Path
from collections.abc import Iterable, Mapping
from typing import Any

from .active_site import ActiveSiteDefinition
from .analyzer import AnalyzerConfig
from .criteria import Criteria, DistanceHysteresisParams
from .slab import HostDefinition
from .sites import SiteDefinition


ConfigInput = str | Path | Mapping[str, Any]

ALLOWED_GEOMETRIC_SITE_SIGNATURE_MODES = {"none", "type", "id", "type+id", "detailed"}

LEGACY_TOP_LEVEL_KEYS = {"slab", "slab_definition", "site", "reactive_site"}
LEGACY_ANALYZER_KEYS = {
    "site_signature_mode",
    "reactive_site_record_states",
    "reactive_site_record_events",
    "reactive_site_record_joint_reactions",
    "reactive_site_streaming",
}


@dataclass(frozen=True)
class RunOverrides:
    traj: str | None = None
    out_prefix: str | None = None
    stride: int | None = None
    max_frames: int | None = None
    progress_with_total: bool | None = None

    geometric_site_file: str | None = None
    geometric_site_index_base: int | None = None
    geometric_site_strict: bool | None = None

    active_site_file: str | None = None
    active_site_index_base: int | None = None
    active_site_strict: bool | None = None


@dataclass(frozen=True)
class PreparedRunConfig:
    traj: str
    out_prefix: str
    stride: int
    max_frames: int | None
    progress_with_total: bool
    criteria: Criteria
    slab_def: HostDefinition
    analyzer_config: AnalyzerConfig


def _load_yaml_file(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Loading YAML config requires PyYAML. Please `pip install pyyaml`."
        ) from e

    with path.open("r", encoding="utf-8") as fh:
        obj = yaml.safe_load(fh) or {}
    if not isinstance(obj, Mapping):
        raise ValueError(f"YAML root must be a mapping: {path}")
    return dict(obj)


def _to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"1", "true", "yes", "y", "on"}:
            return True
        if text in {"0", "false", "no", "n", "off", ""}:
            return False
    return bool(value)


def _to_int_or_none(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    return int(value)


def _get_section(parent: Mapping[str, Any], key: str) -> dict[str, Any]:
    value = parent.get(key, {})
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"Section '{key}' must be a mapping/object.")
    return dict(value)


def _as_list(value: Any, field_name: str) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{field_name} must be a list, got string.")
    if isinstance(value, Mapping):
        raise ValueError(f"{field_name} must be a list, got mapping.")
    if not isinstance(value, Iterable):
        raise ValueError(f"{field_name} must be an iterable.")
    return list(value)


def _resolve_path(raw: str, base_dir: Path | None) -> str:
    path = Path(raw).expanduser()
    if path.is_absolute() or base_dir is None:
        return str(path)
    return str((base_dir / path).resolve())


def _normalize_host_index(raw: Any, *, index_base: int, field_name: str) -> int:
    try:
        value = int(raw)
    except Exception as e:
        raise ValueError(f"{field_name} contains non-integer index: {raw!r}") from e

    if index_base == 0:
        if value < 0:
            raise ValueError(f"{field_name} contains invalid 0-based index {value} (must be >= 0).")
        return value
    if index_base == 1:
        if value < 1:
            raise ValueError(f"{field_name} contains invalid 1-based index {value} (must be >= 1).")
        return value - 1
    raise ValueError(f"{field_name} index_base must be 0 or 1, got {index_base}.")


def _parse_host_indices(raw: Any, *, index_base: int, field_name: str) -> frozenset[int]:
    if raw is None:
        return frozenset()

    tokens: list[Any]
    if isinstance(raw, str):
        text = raw.strip()
        if not text:
            return frozenset()
        if text.startswith("[") and text.endswith("]"):
            text = text[1:-1]
        tokens = [part.strip() for part in text.split(",") if part.strip()]
    else:
        tokens = _as_list(raw, field_name)

    out: set[int] = set()
    for token in tokens:
        if isinstance(token, str):
            part = token.strip()
            if not part:
                continue
            if "-" in part:
                start_raw, stop_raw = [x.strip() for x in part.split("-", 1)]
                if not start_raw or not stop_raw:
                    raise ValueError(f"{field_name} contains invalid range token: {token!r}")
                start = _normalize_host_index(start_raw, index_base=index_base, field_name=field_name)
                stop = _normalize_host_index(stop_raw, index_base=index_base, field_name=field_name)
                if start > stop:
                    raise ValueError(f"{field_name} contains descending range: {token!r}")
                out.update(range(start, stop + 1))
                continue
        out.add(_normalize_host_index(token, index_base=index_base, field_name=field_name))
    return frozenset(sorted(out))


def load_raw_config(config: ConfigInput) -> tuple[dict[str, Any], Path | None]:
    if isinstance(config, Mapping):
        return deepcopy(dict(config)), None

    path = Path(config).expanduser()
    return _load_yaml_file(path), path.resolve().parent


def apply_overrides(cfg: Mapping[str, Any], overrides: RunOverrides | None) -> dict[str, Any]:
    out = deepcopy(dict(cfg))
    if overrides is None:
        return out

    run_cfg = out.setdefault("run", {})
    if not isinstance(run_cfg, dict):
        raise ValueError("Section 'run' must be a mapping/object.")

    if overrides.traj is not None:
        run_cfg["traj"] = overrides.traj
    if overrides.out_prefix is not None:
        run_cfg["out_prefix"] = overrides.out_prefix
    if overrides.stride is not None:
        run_cfg["stride"] = int(overrides.stride)
    if overrides.max_frames is not None:
        run_cfg["max_frames"] = int(overrides.max_frames)
    if overrides.progress_with_total is not None:
        run_cfg["progress_with_total"] = bool(overrides.progress_with_total)

    geometric_site_cfg = out.setdefault("geometric_site", {})
    if not isinstance(geometric_site_cfg, dict):
        raise ValueError("Section 'geometric_site' must be a mapping/object.")
    if overrides.geometric_site_file is not None:
        geometric_site_cfg["enabled"] = True
        geometric_site_cfg["file"] = overrides.geometric_site_file
    if overrides.geometric_site_index_base is not None:
        geometric_site_cfg["index_base"] = int(overrides.geometric_site_index_base)
    if overrides.geometric_site_strict is not None:
        geometric_site_cfg["strict_index_validation"] = bool(overrides.geometric_site_strict)

    active_site_cfg = out.setdefault("active_site", {})
    if not isinstance(active_site_cfg, dict):
        raise ValueError("Section 'active_site' must be a mapping/object.")
    if overrides.active_site_file is not None:
        active_site_cfg["enabled"] = True
        active_site_cfg["file"] = overrides.active_site_file
    if overrides.active_site_index_base is not None:
        active_site_cfg["index_base"] = int(overrides.active_site_index_base)
    if overrides.active_site_strict is not None:
        active_site_cfg["strict_core_validation"] = bool(overrides.active_site_strict)

    return out


def _reject_legacy_keys(cfg: Mapping[str, Any]) -> None:
    legacy_top = sorted(key for key in LEGACY_TOP_LEVEL_KEYS if key in cfg)
    if legacy_top:
        raise ValueError(
            "Legacy top-level config keys are no longer supported: "
            f"{legacy_top}. Please use 'host_definition', 'geometric_site', and 'active_site'."
        )

    analyzer = _get_section(cfg, "analyzer")
    legacy_analyzer = sorted(key for key in LEGACY_ANALYZER_KEYS if key in analyzer)
    if legacy_analyzer:
        raise ValueError(
            "Legacy analyzer keys are no longer supported: "
            f"{legacy_analyzer}. Please use 'geometric_site_signature_mode' and "
            "'active_site_record_*' / 'active_site_streaming'."
        )


def _parse_run_config(cfg: Mapping[str, Any], *, base_dir: Path | None) -> tuple[str, str, int, int | None, bool]:
    run = _get_section(cfg, "run")
    traj_raw = str(run.get("traj", "")).strip()
    if not traj_raw:
        raise ValueError("Missing run.traj.")

    out_prefix = str(run.get("out_prefix", "out")).strip() or "out"
    stride = int(run.get("stride", 1))
    if stride < 1:
        raise ValueError("run.stride must be >= 1.")

    max_frames = _to_int_or_none(run.get("max_frames"))
    if max_frames is not None and max_frames < 0:
        raise ValueError("run.max_frames must be >= 0 or null.")

    progress_with_total = _to_bool(run.get("progress_with_total", True))
    traj = _resolve_path(traj_raw, base_dir)
    return traj, out_prefix, stride, max_frames, progress_with_total


def _parse_criteria(cfg: Mapping[str, Any]) -> Criteria:
    section = _get_section(cfg, "criteria")
    cov = _get_section(section, "cov")
    ads = _get_section(section, "ads")
    slab = _get_section(section, "slab")

    return Criteria(
        persist=int(section.get("persist", 5)),
        cooldown=int(section.get("cooldown", 3)),
        cov=DistanceHysteresisParams(
            float(cov.get("form", 1.15)),
            float(cov.get("break", 1.25)),
            _to_int_or_none(cov.get("max_neighbors", None)),
        ),
        ads=DistanceHysteresisParams(
            float(ads.get("form", 1.25)),
            float(ads.get("break", 1.40)),
            _to_int_or_none(ads.get("max_neighbors", None)),
        ),
        slab=DistanceHysteresisParams(
            float(slab.get("form", 1.35)),
            float(slab.get("break", 1.55)),
            _to_int_or_none(slab.get("max_neighbors", 12)),
        ),
    )


def _parse_host_definition(cfg: Mapping[str, Any]) -> HostDefinition:
    section = _get_section(cfg, "host_definition")
    if not section:
        return HostDefinition()

    mode = str(section.get("mode", "")).strip().lower()
    indices_raw = section.get("indices")
    elements_raw = section.get("elements")
    invert = _to_bool(section.get("invert", False))
    index_base = int(section.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("host_definition.index_base must be 0 or 1.")

    if not mode:
        if indices_raw not in (None, []):
            mode = "indices"
        elif elements_raw not in (None, []):
            mode = "elements"
        else:
            return HostDefinition()

    if mode == "indices":
        indices = _parse_host_indices(
            indices_raw,
            index_base=index_base,
            field_name="host_definition.indices",
        )
        return HostDefinition(indices=indices, invert=invert)

    if mode == "elements":
        elements = frozenset(str(x) for x in _as_list(elements_raw, "host_definition.elements"))
        return HostDefinition(elements=elements, invert=invert)

    raise ValueError("host_definition.mode must be 'elements' or 'indices'.")


def _parse_geometric_site_definition(cfg: Mapping[str, Any], *, base_dir: Path | None) -> SiteDefinition | None:
    section = _get_section(cfg, "geometric_site")
    if not _to_bool(section.get("enabled", False)):
        return None

    site_file = str(section.get("file", "")).strip()
    if not site_file:
        raise ValueError("geometric_site.enabled=true but geometric_site.file is missing.")

    index_base = int(section.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("geometric_site.index_base must be 0 or 1.")

    strict = _to_bool(section.get("strict_index_validation", True))
    return SiteDefinition.from_yaml(
        _resolve_path(site_file, base_dir),
        index_base=index_base,
        strict_index_validation=strict,
    )


def _parse_active_site_definition(cfg: Mapping[str, Any], *, base_dir: Path | None) -> ActiveSiteDefinition | None:
    section = _get_section(cfg, "active_site")
    if not _to_bool(section.get("enabled", False)):
        return None

    site_file = str(section.get("file", "")).strip()
    if not site_file:
        raise ValueError("active_site.enabled=true but active_site.file is missing.")

    index_base = int(section.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("active_site.index_base must be 0 or 1.")

    strict = _to_bool(section.get("strict_core_validation", True))
    return ActiveSiteDefinition.from_yaml(
        _resolve_path(site_file, base_dir),
        index_base=index_base,
        strict_core_validation=strict,
    )


def _parse_analyzer_config(
    cfg: Mapping[str, Any],
    *,
    geometric_site_definition: SiteDefinition | None,
    active_site_definition: ActiveSiteDefinition | None,
    out_prefix: str,
    strict_fields: bool = False,
) -> AnalyzerConfig:
    analyzer = _get_section(cfg, "analyzer")
    translated: dict[str, Any] = dict(analyzer)

    mode_raw = translated.get("geometric_site_signature_mode", "none")
    mode = str(mode_raw).strip().lower()
    if mode not in ALLOWED_GEOMETRIC_SITE_SIGNATURE_MODES:
        raise ValueError(
            f"Invalid analyzer.geometric_site_signature_mode={mode_raw!r}; "
            f"allowed: {sorted(ALLOWED_GEOMETRIC_SITE_SIGNATURE_MODES)}"
        )

    if geometric_site_definition is None:
        translated["geometric_site_signature_mode"] = "none"
    else:
        translated["geometric_site_signature_mode"] = mode

    translated["geometric_site_definition"] = geometric_site_definition
    translated["active_site_definition"] = active_site_definition

    if _to_bool(translated.get("record_frame_species", False)) and _to_bool(
        translated.get("frame_species_streaming", True)
    ):
        if not translated.get("frame_species_path"):
            translated["frame_species_path"] = f"{out_prefix}_frame_species.csv"

    valid_fields = {field.name for field in fields(AnalyzerConfig)}
    unknown = sorted(key for key in translated if key not in valid_fields)
    if unknown and strict_fields:
        raise ValueError(f"Unknown analyzer keys: {unknown}")

    filtered = {key: value for key, value in translated.items() if key in valid_fields}
    return AnalyzerConfig(**filtered)


def build_prepared_config(
    config: ConfigInput,
    *,
    overrides: RunOverrides | None = None,
    base_dir: str | Path | None = None,
    strict_analyzer_fields: bool = False,
) -> PreparedRunConfig:
    cfg, inferred_base = load_raw_config(config)
    cfg = apply_overrides(cfg, overrides)
    _reject_legacy_keys(cfg)

    resolved_base = inferred_base if base_dir is None else Path(base_dir).expanduser().resolve()

    traj, out_prefix, stride, max_frames, progress_with_total = _parse_run_config(cfg, base_dir=resolved_base)
    criteria = _parse_criteria(cfg)
    host_definition = _parse_host_definition(cfg)
    geometric_site_definition = _parse_geometric_site_definition(cfg, base_dir=resolved_base)
    active_site_definition = _parse_active_site_definition(cfg, base_dir=resolved_base)
    analyzer_config = _parse_analyzer_config(
        cfg,
        geometric_site_definition=geometric_site_definition,
        active_site_definition=active_site_definition,
        out_prefix=out_prefix,
        strict_fields=strict_analyzer_fields,
    )

    return PreparedRunConfig(
        traj=traj,
        out_prefix=out_prefix,
        stride=stride,
        max_frames=max_frames,
        progress_with_total=progress_with_total,
        criteria=criteria,
        slab_def=host_definition,
        analyzer_config=analyzer_config,
    )
