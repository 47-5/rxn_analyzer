"""YAML-driven graph postprocess entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx

from .context import run_context_mode
from .focus import run_focus_mode
from .story import run_story_mode


def load_postprocess_config(path: str) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError("Loading graph postprocess YAML requires PyYAML.") from e

    with open(path, "r", encoding="utf-8") as f:
        obj = yaml.safe_load(f) or {}
    if not isinstance(obj, dict):
        raise ValueError("Postprocess config root must be a YAML mapping.")
    return obj


def run_postprocess_from_config(config_path: str) -> nx.Graph:
    cfg = load_postprocess_config(config_path)
    mode = str(cfg.get("mode", "focus") or "focus").strip().lower()

    input_cfg = cfg.get("input", {})
    output_cfg = cfg.get("output", {})
    graph_path = str(input_cfg.get("graph", "")).strip()
    if not graph_path:
        raise ValueError("Postprocess config must define input.graph.")

    out_graph_path = str(output_cfg.get("graph", "")).strip()
    if not out_graph_path:
        raise ValueError("Postprocess config must define output.graph.")

    score_csv_path = str(output_cfg.get("scores_csv", "")).strip() or None

    print(f"[read] {graph_path}")
    graph = nx.read_graphml(graph_path)

    if mode == "focus":
        result = run_focus_mode(graph, cfg, score_csv_path=score_csv_path)
    elif mode == "context":
        result = run_context_mode(graph, cfg)
    elif mode == "story":
        result = run_story_mode(graph, cfg)
    else:
        raise ValueError(f"Unknown postprocess mode: {mode}")

    out_dir = Path(out_graph_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[write] {out_graph_path}")
    nx.write_graphml(result, out_graph_path)
    return result


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser("graph postprocess")
    parser.add_argument("--config", required=True, help="YAML postprocess config")
    args = parser.parse_args()
    run_postprocess_from_config(args.config)


if __name__ == "__main__":
    main()
