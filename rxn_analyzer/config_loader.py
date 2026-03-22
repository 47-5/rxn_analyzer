from __future__ import annotations

import warnings
from copy import deepcopy
from dataclasses import dataclass, fields
from pathlib import Path
from collections.abc import Iterable, Mapping
from typing import Any

from .analyzer import AnalyzerConfig
from .criteria import Criteria, DistanceHysteresisParams
from .slab import SlabDefinition
from .sites import SiteDefinition
from .active_site import ActiveSiteDefinition


ConfigInput = str | Path | Mapping[str, Any]

_ALLOWED_SITE_SIGNATURE_MODES = {"none", "type", "id", "type+id", "detailed"}


def _get_section_alias(parent: Mapping[str, Any], *keys: str) -> dict[str, Any]:
    for key in keys:
        if key in parent:
            return _get_section(parent, key)
    return {}


@dataclass(frozen=True)
class RunOverrides:
    traj: str | None = None
    out_prefix: str | None = None
    stride: int | None = None
    max_frames: int | None = None
    progress_with_total: bool | None = None

    site_file: str | None = None
    site_index_base: int | None = None
    site_strict: bool | None = None
    geometric_site_file: str | None = None
    geometric_site_index_base: int | None = None
    geometric_site_strict: bool | None = None
    reactive_site_file: str | None = None
    reactive_site_index_base: int | None = None
    reactive_site_strict: bool | None = None
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
    slab_def: SlabDefinition
    analyzer_config: AnalyzerConfig


def _load_yaml_file(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Loading YAML config requires PyYAML. Please `pip install pyyaml`."
        ) from e

    with path.open("r", encoding="utf-8") as f:
        obj = yaml.safe_load(f) or {}

    if not isinstance(obj, Mapping):
        raise ValueError(f"YAML root must be a mapping: {path}")
    return dict(obj)


def _to_bool(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(v)
    if isinstance(v, str):
        s = v.strip().lower()
        if s in {"1", "true", "yes", "y", "on"}:
            return True
        if s in {"0", "false", "no", "n", "off", ""}:
            return False
    return bool(v)


def _to_int_or_none(v: Any) -> int | None:
    if v is None:
        return None
    if isinstance(v, str) and not v.strip():
        return None
    return int(v)


def _get_section(parent: Mapping[str, Any], key: str) -> dict[str, Any]:
    v = parent.get(key, {})
    if v is None:
        return {}
    if not isinstance(v, Mapping):
        raise ValueError(f"Section '{key}' must be a mapping/object.")
    return dict(v)


def _as_list(v: Any, field_name: str) -> list[Any]:
    if v is None:
        return []
    if isinstance(v, (str, bytes)):
        raise ValueError(f"{field_name} must be a list, got string.")
    if isinstance(v, Mapping):
        raise ValueError(f"{field_name} must be a list, got mapping.")
    if not isinstance(v, Iterable):
        raise ValueError(f"{field_name} must be an iterable.")
    return list(v)


def _normalize_host_index(raw: Any, *, index_base: int, field_name: str) -> int:
    try:
        v = int(raw)
    except Exception as e:
        raise ValueError(f"{field_name} contains non-integer index: {raw!r}") from e

    if index_base == 0:
        if v < 0:
            raise ValueError(f"{field_name} contains invalid 0-based index {v} (must be >= 0).")
        return v
    if index_base == 1:
        if v < 1:
            raise ValueError(f"{field_name} contains invalid 1-based index {v} (must be >= 1).")
        return v - 1
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


def _resolve_path(raw: str, base_dir: Path | None) -> str:
    p = Path(raw).expanduser()
    if p.is_absolute() or base_dir is None:
        return str(p)
    return str((base_dir / p).resolve())


def load_raw_config(config: ConfigInput) -> tuple[dict[str, Any], Path | None]:
    """
    Returns:
      - cfg dict
      - inferred base_dir (if config is path), else None
    """
    if isinstance(config, Mapping):
        return deepcopy(dict(config)), None

    path = Path(config).expanduser()
    cfg = _load_yaml_file(path)
    return cfg, path.resolve().parent


def apply_overrides(cfg: Mapping[str, Any], overrides: RunOverrides | None) -> dict[str, Any]:
    out = deepcopy(dict(cfg))
    if overrides is None:
        return out

    run_cfg = out.get("run")
    if not isinstance(run_cfg, dict):
        run_cfg = {}
        out["run"] = run_cfg

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

    site_cfg = out.get("geometric_site")
    if not isinstance(site_cfg, dict):
        site_cfg = out.get("site")
    if not isinstance(site_cfg, dict):
        site_cfg = {}
    out["site"] = site_cfg
    out["geometric_site"] = site_cfg

    effective_site_file = overrides.geometric_site_file if overrides.geometric_site_file is not None else overrides.site_file
    effective_site_index_base = (
        overrides.geometric_site_index_base if overrides.geometric_site_index_base is not None else overrides.site_index_base
    )
    effective_site_strict = (
        overrides.geometric_site_strict if overrides.geometric_site_strict is not None else overrides.site_strict
    )

    if effective_site_file is not None:
        site_cfg["enabled"] = True
        site_cfg["file"] = effective_site_file
    if effective_site_index_base is not None:
        site_cfg["index_base"] = int(effective_site_index_base)
    if effective_site_strict is not None:
        site_cfg["strict_index_validation"] = bool(effective_site_strict)

    reactive_site_cfg = out.get("active_site")
    if not isinstance(reactive_site_cfg, dict):
        reactive_site_cfg = out.get("reactive_site")
    if not isinstance(reactive_site_cfg, dict):
        reactive_site_cfg = {}
    out["reactive_site"] = reactive_site_cfg
    out["active_site"] = reactive_site_cfg

    effective_reactive_file = (
        overrides.active_site_file if overrides.active_site_file is not None else overrides.reactive_site_file
    )
    effective_reactive_index_base = (
        overrides.active_site_index_base
        if overrides.active_site_index_base is not None
        else overrides.reactive_site_index_base
    )
    effective_reactive_strict = (
        overrides.active_site_strict if overrides.active_site_strict is not None else overrides.reactive_site_strict
    )

    if effective_reactive_file is not None:
        reactive_site_cfg["enabled"] = True
        reactive_site_cfg["file"] = effective_reactive_file
    if effective_reactive_index_base is not None:
        reactive_site_cfg["index_base"] = int(effective_reactive_index_base)
    if effective_reactive_strict is not None:
        reactive_site_cfg["strict_core_validation"] = bool(effective_reactive_strict)

    return out


def enforce_site_consistency(cfg: Mapping[str, Any], *, strict: bool = False) -> dict[str, Any]:
    """
    校验并规范 site 与 analyzer.site_signature_mode 的一致性。

    规则：
    1) site.enabled=false 时：
       - strict=True：若 mode != 'none' 则报错
       - strict=False：自动将 mode 改成 'none'
       - 且强制忽略 site.file，避免误加载
    2) site.enabled=true 时：site.file 必须提供
    3) analyzer.site_signature_mode 必须在允许值内
    """
    out = deepcopy(dict(cfg))

    site_raw = out.get("geometric_site", out.get("site", {}))
    if site_raw is None:
        site_raw = {}
    if not isinstance(site_raw, Mapping):
        raise ValueError("Section 'site' must be a mapping/object.")
    site = dict(site_raw)
    out["site"] = site
    out["geometric_site"] = site

    analyzer_raw = out.get("analyzer", {})
    if analyzer_raw is None:
        analyzer_raw = {}
    if not isinstance(analyzer_raw, Mapping):
        raise ValueError("Section 'analyzer' must be a mapping/object.")
    analyzer = dict(analyzer_raw)
    out["analyzer"] = analyzer

    enabled = _to_bool(site.get("enabled", False))

    if "geometric_site_signature_mode" in analyzer and "site_signature_mode" not in analyzer:
        analyzer["site_signature_mode"] = analyzer["geometric_site_signature_mode"]
    mode_raw = analyzer.get("site_signature_mode", "none")
    mode = str(mode_raw).strip().lower() if mode_raw is not None else "none"
    if not mode:
        mode = "none"

    if mode not in _ALLOWED_SITE_SIGNATURE_MODES:
        raise ValueError(
            f"Invalid analyzer.site_signature_mode={mode_raw!r}, "
            f"allowed: {sorted(_ALLOWED_SITE_SIGNATURE_MODES)}"
        )

    if not enabled:
        if mode != "none":
            msg = (
                "Config conflict: site.enabled=false but "
                f"analyzer.site_signature_mode={mode_raw!r} (must be 'none' when site is disabled)."
            )
            if strict:
                raise ValueError(msg)
            warnings.warn(msg + " Auto-set to 'none'.", RuntimeWarning)
            analyzer["site_signature_mode"] = "none"

        # 关键：禁用时强制忽略 file，防止后续误加载 site
        if site.get("file") not in (None, ""):
            warnings.warn("site.enabled=false; site.file will be ignored.", RuntimeWarning)
        site["file"] = None
        site["strict_index_validation"] = False

    else:
        site_file = site.get("file", None)
        if site_file is None or (isinstance(site_file, str) and not site_file.strip()):
            raise ValueError("site.enabled=true but site.file is missing.")
        if mode == "none":
            warnings.warn(
                "site.enabled=true but analyzer.site_signature_mode='none'; "
                "site will be loaded/validated but not appended to species labels.",
                RuntimeWarning,
            )

    return out


def _build_criteria(cfg: Mapping[str, Any]) -> Criteria:
    c = _get_section(cfg, "criteria")
    cov = _get_section(c, "cov")
    ads = _get_section(c, "ads")
    slab = _get_section(c, "slab")

    return Criteria(
        persist=int(c.get("persist", 5)),
        cooldown=int(c.get("cooldown", 3)),
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


def _build_slab_def(cfg: Mapping[str, Any]) -> SlabDefinition:
    # 支持 slab_definition（推荐）或 slab（兼容之前示例）
    s = _get_section_alias(cfg, "host_definition", "slab_definition")
    if not s:
        s = _get_section_alias(cfg, "host", "slab")

    if not s:
        return SlabDefinition()

    has_def_keys = any(k in s for k in ("mode", "elements", "indices", "invert", "index_base"))
    if not has_def_keys:
        # 如果 top-level slab 没有定义这些字段，视为未配置 slab 定义
        return SlabDefinition()

    mode = str(s.get("mode", "")).strip().lower()
    indices_raw = s.get("indices", None)
    elements_raw = s.get("elements", None)
    invert = _to_bool(s.get("invert", False))
    index_base = int(s.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("host_definition.index_base must be 0 or 1.")

    if not mode:
        if indices_raw not in (None, []):
            mode = "indices"
        elif elements_raw not in (None, []):
            mode = "elements"
        else:
            return SlabDefinition()

    if mode == "indices":
        idx = _parse_host_indices(indices_raw, index_base=index_base, field_name="host_definition.indices")
        return SlabDefinition(indices=idx, invert=invert)

    if mode == "elements":
        elems = [str(x) for x in _as_list(elements_raw, "slab.elements")]
        return SlabDefinition(elements=frozenset(elems), invert=invert)

    raise ValueError(f"Unknown slab mode: {mode}. Use 'elements' or 'indices'.")


def _build_site_def(cfg: Mapping[str, Any], base_dir: Path | None) -> SiteDefinition | None:
    sc = _get_section_alias(cfg, "geometric_site", "site")

    # 关键修复：enabled=false 时，无论 file 是否填写都不加载 site
    enabled = _to_bool(sc.get("enabled", False))
    if not enabled:
        return None

    site_file = sc.get("file", None)
    if not site_file:
        raise ValueError("site.enabled=true but site.file is missing.")

    index_base = int(sc.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("site.index_base must be 0 or 1.")

    strict = _to_bool(sc.get("strict_index_validation", True))
    site_path = _resolve_path(str(site_file), base_dir)

    return SiteDefinition.from_yaml(
        site_path,
        index_base=index_base,
        strict_index_validation=strict,
    )


def _build_reactive_site_def(cfg: Mapping[str, Any], base_dir: Path | None) -> ActiveSiteDefinition | None:
    sc = _get_section_alias(cfg, "active_site", "reactive_site")
    enabled = _to_bool(sc.get("enabled", False))
    if not enabled:
        return None

    site_file = sc.get("file", None)
    if not site_file:
        raise ValueError("reactive_site.enabled=true but reactive_site.file is missing.")

    index_base = int(sc.get("index_base", 0))
    if index_base not in (0, 1):
        raise ValueError("reactive_site.index_base must be 0 or 1.")

    strict = _to_bool(sc.get("strict_core_validation", True))
    site_path = _resolve_path(str(site_file), base_dir)
    return ActiveSiteDefinition.from_yaml(
        site_path,
        index_base=index_base,
        strict_core_validation=strict,
    )


def _build_analyzer_config(
    cfg: Mapping[str, Any],
    *,
    site_def: SiteDefinition | None,
    reactive_site_def: ActiveSiteDefinition | None,
    out_prefix: str,
    strict_analyzer_fields: bool = False,
) -> AnalyzerConfig:
    a = _get_section(cfg, "analyzer")

    # 兼容旧 CLI：allow_suspicious => allow charged + dot, 且不 fallback
    allow_suspicious = _to_bool(a.pop("smiles_allow_suspicious", False))
    if allow_suspicious:
        a["smiles_allow_charged"] = True
        a["smiles_allow_dot"] = True
        a["smiles_fallback_to_formula_if_suspicious"] = False

    if site_def is not None:
        a["site_definition"] = site_def
    if reactive_site_def is not None:
        a["active_site_definition"] = reactive_site_def
    if "geometric_site_signature_mode" in a and "site_signature_mode" not in a:
        a["site_signature_mode"] = a["geometric_site_signature_mode"]
    if "active_site_record_states" in a and "reactive_site_record_states" not in a:
        a["reactive_site_record_states"] = a["active_site_record_states"]
    if "active_site_record_events" in a and "reactive_site_record_events" not in a:
        a["reactive_site_record_events"] = a["active_site_record_events"]
    if "active_site_record_joint_reactions" in a and "reactive_site_record_joint_reactions" not in a:
        a["reactive_site_record_joint_reactions"] = a["active_site_record_joint_reactions"]
    if "active_site_streaming" in a and "reactive_site_streaming" not in a:
        a["reactive_site_streaming"] = a["active_site_streaming"]
    if "reactive_site_streaming" in a and "active_site_streaming" not in a:
        a["active_site_streaming"] = a["reactive_site_streaming"]

    # 对齐旧 CLI：记录 frame species + streaming 且未给 path 时，自动补默认路径
    if _to_bool(a.get("record_frame_species", False)) and _to_bool(a.get("frame_species_streaming", True)):
        if not a.get("frame_species_path"):
            a["frame_species_path"] = f"{out_prefix}_frame_species.csv"

    valid = {f.name for f in fields(AnalyzerConfig)}
    unknown = sorted(k for k in a.keys() if k not in valid)
    if unknown and strict_analyzer_fields:
        raise ValueError(f"Unknown analyzer keys: {unknown}")

    filtered = {k: v for k, v in a.items() if k in valid}
    return AnalyzerConfig(**filtered)


def build_prepared_config(
    config: ConfigInput,
    *,
    overrides: RunOverrides | None = None,
    base_dir: str | Path | None = None,
    strict_analyzer_fields: bool = False,
    strict_site_consistency: bool = False,
) -> PreparedRunConfig:
    """
    只做配置加载和对象构建，不执行分析。
    """
    cfg, inferred_base = load_raw_config(config)
    merged = apply_overrides(cfg, overrides)
    merged = enforce_site_consistency(merged, strict=strict_site_consistency)

    resolved_base: Path | None
    if base_dir is None:
        resolved_base = inferred_base
    else:
        resolved_base = Path(base_dir).expanduser().resolve()

    run = _get_section(merged, "run")
    traj_raw = str(run.get("traj", "")).strip()
    if not traj_raw:
        raise ValueError("Missing run.traj (or override traj).")

    traj = _resolve_path(traj_raw, resolved_base)

    out_prefix = str(run.get("out_prefix", "out"))
    stride = int(run.get("stride", 1))
    if stride < 1:
        raise ValueError("run.stride must be >= 1.")

    max_frames_raw = run.get("max_frames", None)
    max_frames = None if max_frames_raw is None else int(max_frames_raw)
    if max_frames is not None and max_frames < 0:
        raise ValueError("run.max_frames must be >= 0 or null.")
    progress_with_total = _to_bool(run.get("progress_with_total", True))

    criteria = _build_criteria(merged)
    slab_def = _build_slab_def(merged)
    site_def = _build_site_def(merged, resolved_base)
    reactive_site_def = _build_reactive_site_def(merged, resolved_base)
    analyzer_cfg = _build_analyzer_config(
        merged,
        site_def=site_def,
        reactive_site_def=reactive_site_def,
        out_prefix=out_prefix,
        strict_analyzer_fields=strict_analyzer_fields,
    )

    return PreparedRunConfig(
        traj=traj,
        out_prefix=out_prefix,
        stride=stride,
        max_frames=max_frames,
        progress_with_total=progress_with_total,
        criteria=criteria,
        slab_def=slab_def,
        analyzer_config=analyzer_cfg,
    )
