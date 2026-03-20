from __future__ import annotations

import argparse
from collections.abc import Sequence

from tqdm import tqdm

from . import network as net
from .analyzer import ReactionAnalyzer
from .config_loader import (
    ConfigInput,
    PreparedRunConfig,
    RunOverrides,
    build_prepared_config,
)
from .io import frames


def _count_frames(traj: str, stride: int) -> int:
    return sum(1 for _ in frames(traj, stride=stride))


def execute_prepared_config(
    prepared: PreparedRunConfig,
    *,
    show_progress: bool = True,
    progress_with_total: bool = False,
    reset_node_mapping: bool = True,
) -> ReactionAnalyzer:
    if reset_node_mapping:
        net.reset_node_id_mapping()

    analyzer = ReactionAnalyzer(
        criteria=prepared.criteria,
        slab_def=prepared.slab_def,
        config=prepared.analyzer_config,
    )

    stream = frames(prepared.traj, stride=prepared.stride)
    if show_progress:
        if progress_with_total:
            total = _count_frames(prepared.traj, prepared.stride)
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
    site_file: str | None = None,
    site_index_base: int | None = None,
    site_strict: bool | None = None,
    reactive_site_file: str | None = None,
    reactive_site_index_base: int | None = None,
    reactive_site_strict: bool | None = None,
    base_dir: str | None = None,
    show_progress: bool = True,
    progress_with_total: bool = False,
    reset_node_mapping: bool = True,
    strict_analyzer_fields: bool = False,
    strict_site_consistency: bool = False,
) -> ReactionAnalyzer:
    overrides = RunOverrides(
        traj=traj,
        out_prefix=out_prefix,
        stride=stride,
        max_frames=max_frames,
        site_file=site_file,
        site_index_base=site_index_base,
        site_strict=site_strict,
        reactive_site_file=reactive_site_file,
        reactive_site_index_base=reactive_site_index_base,
        reactive_site_strict=reactive_site_strict,
    )

    prepared = build_prepared_config(
        config,
        overrides=overrides,
        base_dir=base_dir,
        strict_analyzer_fields=strict_analyzer_fields,
        strict_site_consistency=strict_site_consistency,
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
    p.add_argument("--site-file", default=None, help="Override site.file")
    p.add_argument("--site-index-base", type=int, choices=[0, 1], default=None, help="Override site.index_base")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--site-strict", dest="site_strict", action="store_true", help="Override strict_index_validation=true")
    g.add_argument("--no-site-strict", dest="site_strict", action="store_false", help="Override strict_index_validation=false")
    p.set_defaults(site_strict=None)

    p.add_argument("--reactive-site-file", default=None, help="Override reactive_site.file")
    p.add_argument(
        "--reactive-site-index-base",
        type=int,
        choices=[0, 1],
        default=None,
        help="Override reactive_site.index_base",
    )
    g2 = p.add_mutually_exclusive_group()
    g2.add_argument(
        "--reactive-site-strict",
        dest="reactive_site_strict",
        action="store_true",
        help="Override reactive_site.strict_core_validation=true",
    )
    g2.add_argument(
        "--no-reactive-site-strict",
        dest="reactive_site_strict",
        action="store_false",
        help="Override reactive_site.strict_core_validation=false",
    )
    p.set_defaults(reactive_site_strict=None)

    p.add_argument(
        "--base-dir",
        default=None,
        help="Base dir for relative paths (default: directory of --config).",
    )

    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar")
    p.add_argument("--progress-with-total", action="store_true", help="Pre-scan trajectory to show total progress")
    p.add_argument("--no-reset-node-mapping", action="store_true", help="Do not call network.reset_node_id_mapping()")

    p.add_argument(
        "--strict-analyzer-fields",
        action="store_true",
        help="Raise error if analyzer section contains unknown keys.",
    )
    p.add_argument(
        "--strict-site-consistency",
        action="store_true",
        help="Raise error on site/analyzer consistency conflicts (instead of auto-fix).",
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
        site_file=args.site_file,
        site_index_base=args.site_index_base,
        site_strict=args.site_strict,
        reactive_site_file=args.reactive_site_file,
        reactive_site_index_base=args.reactive_site_index_base,
        reactive_site_strict=args.reactive_site_strict,
        base_dir=args.base_dir,
        show_progress=(not args.no_progress),
        progress_with_total=bool(args.progress_with_total),
        reset_node_mapping=(not args.no_reset_node_mapping),
        strict_analyzer_fields=bool(args.strict_analyzer_fields),
        strict_site_consistency=bool(args.strict_site_consistency),
    )
    print(f"[done] outputs written with prefix: {analyzer.config.frame_species_path or args.out_prefix or 'out'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
