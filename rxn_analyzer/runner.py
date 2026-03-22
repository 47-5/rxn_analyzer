from __future__ import annotations

import argparse
from collections.abc import Sequence

from tqdm import tqdm

from . import graph as net
from .analyzer import ReactionAnalyzer
from .config_loader import (
    ConfigInput,
    PreparedRunConfig,
    RunOverrides,
    build_prepared_config,
)
from .io import frames


def _count_frames(traj: str, stride: int, max_frames: int | None = None) -> int:
    total = sum(1 for _ in frames(traj, stride=stride))
    if max_frames is not None:
        return min(total, int(max_frames))
    return total


def execute_prepared_config(
    prepared: PreparedRunConfig,
    *,
    show_progress: bool = True,
    progress_with_total: bool | None = None,
    reset_node_mapping: bool = True,
) -> ReactionAnalyzer:
    if reset_node_mapping:
        net.reset_node_id_mapping()

    analyzer = ReactionAnalyzer(
        criteria=prepared.criteria,
        slab_def=prepared.slab_def,
        config=prepared.analyzer_config,
        out_prefix=prepared.out_prefix,
    )

    use_progress_with_total = (
        prepared.progress_with_total if progress_with_total is None else bool(progress_with_total)
    )

    stream = frames(prepared.traj, stride=prepared.stride)
    if show_progress:
        if use_progress_with_total:
            total = _count_frames(prepared.traj, prepared.stride, prepared.max_frames)
            stream = tqdm(stream, total=total, desc="Analyzing", unit="frame")
        else:
            stream = tqdm(stream, desc="Analyzing", unit="frame")

    frame_count = 0
    for _k, atoms in stream:
        if prepared.max_frames is not None and frame_count >= prepared.max_frames:
            break
        analyzer.process_frame(frame_count, atoms)
        frame_count += 1

    analyzer.write_outputs(prepared.out_prefix)
    return analyzer


def run_from_yaml(
    config: ConfigInput,
    *,
    traj: str | None = None,
    out_prefix: str | None = None,
    stride: int | None = None,
    max_frames: int | None = None,
    geometric_site_file: str | None = None,
    geometric_site_index_base: int | None = None,
    geometric_site_strict: bool | None = None,
    active_site_file: str | None = None,
    active_site_index_base: int | None = None,
    active_site_strict: bool | None = None,
    base_dir: str | None = None,
    show_progress: bool = True,
    progress_with_total: bool | None = None,
    reset_node_mapping: bool = True,
    strict_analyzer_fields: bool = False,
) -> ReactionAnalyzer:
    overrides = RunOverrides(
        traj=traj,
        out_prefix=out_prefix,
        stride=stride,
        max_frames=max_frames,
        progress_with_total=progress_with_total,
        geometric_site_file=geometric_site_file,
        geometric_site_index_base=geometric_site_index_base,
        geometric_site_strict=geometric_site_strict,
        active_site_file=active_site_file,
        active_site_index_base=active_site_index_base,
        active_site_strict=active_site_strict,
    )

    prepared = build_prepared_config(
        config,
        overrides=overrides,
        base_dir=base_dir,
        strict_analyzer_fields=strict_analyzer_fields,
    )

    return execute_prepared_config(
        prepared,
        show_progress=show_progress,
        progress_with_total=progress_with_total,
        reset_node_mapping=reset_node_mapping,
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser("Run rxn_analyzer from YAML config")

    p.add_argument("--config", required=True, help="Path to analyzer_config.yaml")

    # run overrides
    p.add_argument("--traj", default=None, help="Override run.traj")
    p.add_argument("--out-prefix", default=None, help="Override run.out_prefix")
    p.add_argument("--stride", type=int, default=None, help="Override run.stride")
    p.add_argument("--max-frames", type=int, default=None, help="Override run.max_frames")

    # site overrides
    p.add_argument("--geometric-site-file", default=None, help="Override geometric_site.file")
    p.add_argument("--geometric-site-index-base", type=int, choices=[0, 1], default=None, help="Override geometric_site.index_base")
    g_geo = p.add_mutually_exclusive_group()
    g_geo.add_argument("--geometric-site-strict", dest="geometric_site_strict", action="store_true", help="Override geometric_site.strict_index_validation=true")
    g_geo.add_argument("--no-geometric-site-strict", dest="geometric_site_strict", action="store_false", help="Override geometric_site.strict_index_validation=false")
    p.set_defaults(geometric_site_strict=None)

    p.add_argument("--active-site-file", default=None, help="Override active_site.file")
    p.add_argument(
        "--active-site-index-base",
        type=int,
        choices=[0, 1],
        default=None,
        help="Override active_site.index_base",
    )
    g2a = p.add_mutually_exclusive_group()
    g2a.add_argument(
        "--active-site-strict",
        dest="active_site_strict",
        action="store_true",
        help="Override active_site.strict_core_validation=true",
    )
    g2a.add_argument(
        "--no-active-site-strict",
        dest="active_site_strict",
        action="store_false",
        help="Override active_site.strict_core_validation=false",
    )
    p.set_defaults(active_site_strict=None)

    p.add_argument(
        "--base-dir",
        default=None,
        help="Base dir for relative paths (default: directory of --config).",
    )

    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar")
    g3 = p.add_mutually_exclusive_group()
    g3.add_argument(
        "--progress-with-total",
        dest="progress_with_total",
        action="store_true",
        help="Override run.progress_with_total=true.",
    )
    g3.add_argument(
        "--no-progress-total",
        dest="progress_with_total",
        action="store_false",
        help="Override run.progress_with_total=false.",
    )
    p.set_defaults(progress_with_total=None)
    p.add_argument("--no-reset-node-mapping", action="store_true", help="Do not call network.reset_node_id_mapping()")

    p.add_argument(
        "--strict-analyzer-fields",
        action="store_true",
        help="Raise error if analyzer section contains unknown keys.",
    )

    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    analyzer = run_from_yaml(
        args.config,
        traj=args.traj,
        out_prefix=args.out_prefix,
        stride=args.stride,
        max_frames=args.max_frames,
        geometric_site_file=args.geometric_site_file,
        geometric_site_index_base=args.geometric_site_index_base,
        geometric_site_strict=args.geometric_site_strict,
        active_site_file=args.active_site_file,
        active_site_index_base=args.active_site_index_base,
        active_site_strict=args.active_site_strict,
        base_dir=args.base_dir,
        show_progress=(not args.no_progress),
        progress_with_total=args.progress_with_total,
        reset_node_mapping=(not args.no_reset_node_mapping),
        strict_analyzer_fields=bool(args.strict_analyzer_fields),
    )
    print(f"[done] outputs written with prefix: {analyzer.config.frame_species_path or args.out_prefix or 'out'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
