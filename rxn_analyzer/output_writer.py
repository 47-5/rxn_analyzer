from __future__ import annotations
import csv
import json

import networkx as nx

from .network import (
    ensure_graph,
    add_transform,
    summarize_reversible_reactions,
    finalize_graph_for_export,
)


class OutputWriter:
    def write_all(
        self,
        out_prefix: str,
        transform_events: list,
        bond_events: list,
        graph: nx.DiGraph,
        frame_species_logger,
        reaction_summary_include_frames: bool = True,
    ) -> None:
        tf_path = f"{out_prefix}_events_transform.csv"
        print(f"[write] {tf_path}")
        with open(tf_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=["event_id", "frame", "type", "reactants", "products", "delta_counts", "evidence_bonds", "confidence"],
            )
            w.writeheader()
            for e in transform_events:
                w.writerow(
                    {
                        "event_id": e.event_id,
                        "frame": e.frame,
                        "type": e.type,
                        "reactants": json.dumps(e.reactants, ensure_ascii=False),
                        "products": json.dumps(e.products, ensure_ascii=False),
                        "delta_counts": json.dumps(e.delta_counts, ensure_ascii=False),
                        "evidence_bonds": json.dumps(e.evidence_bonds, ensure_ascii=False),
                        "confidence": e.confidence,
                    }
                )

        be_path = f"{out_prefix}_events_bonds.csv"
        print(f"[write] {be_path}")
        with open(be_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=["event_id", "frame", "edge_type", "action", "i", "j", "distance", "threshold_form", "threshold_break"],
            )
            w.writeheader()
            for e in bond_events:
                w.writerow(
                    {
                        "event_id": e.event_id,
                        "frame": e.frame,
                        "edge_type": e.edge.type.value,
                        "action": e.action,
                        "i": e.edge.i,
                        "j": e.edge.j,
                        "distance": e.evidence.distance,
                        "threshold_form": e.evidence.threshold_form,
                        "threshold_break": e.evidence.threshold_break,
                    }
                )

        print(f"[write] {out_prefix}_network.graphml")
        finalize_graph_for_export(graph)
        nx.write_graphml(graph, f"{out_prefix}_network.graphml")

        legacy = ensure_graph()
        for e in transform_events:
            # 兼容你已修改/未修改 network.add_transform 的两种签名
            try:
                add_transform(legacy, e.reactants, e.products, e.type, frame=e.frame)
            except TypeError:
                add_transform(legacy, e.reactants, e.products, e.type)

        # 兼容你已修改/未修改 network.summarize_reversible_reactions 的两种签名
        try:
            lines = summarize_reversible_reactions(
                legacy,
                use="orig_id",
                min_total_weight=1,
                sort_by="total",
                include_frames=bool(reaction_summary_include_frames),
                unique_frames=False,
            )
        except TypeError:
            lines = summarize_reversible_reactions(
                legacy,
                use="orig_id",
                min_total_weight=1,
                sort_by="total",
            )

        summary_path = f"{out_prefix}_reactions_summary.tsv"
        print(f"[write] {summary_path}")
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

        if frame_species_logger is not None:
            frame_species_logger.finalize(out_prefix)
