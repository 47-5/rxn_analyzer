from __future__ import annotations

import csv
import json

import networkx as nx

from .graph import (
    add_transform,
    ensure_graph,
    finalize_graph_for_export,
    summarize_reversible_reactions,
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
        self._write_transform_events(out_prefix, transform_events)
        self._write_bond_events(out_prefix, bond_events)
        self._write_main_graph(out_prefix, graph)
        self._write_reaction_summary(
            out_prefix,
            transform_events,
            include_frames=bool(reaction_summary_include_frames),
        )
        if frame_species_logger is not None:
            frame_species_logger.finalize(out_prefix)

    def _write_transform_events(self, out_prefix: str, transform_events: list) -> None:
        path = f"{out_prefix}_events_transform.csv"
        print(f"[write] {path}")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "event_id",
                    "frame",
                    "type",
                    "reactants",
                    "products",
                    "delta_counts",
                    "evidence_bonds",
                    "confidence",
                ],
            )
            writer.writeheader()
            for event in transform_events:
                writer.writerow(
                    {
                        "event_id": event.event_id,
                        "frame": event.frame,
                        "type": event.type,
                        "reactants": json.dumps(event.reactants, ensure_ascii=False),
                        "products": json.dumps(event.products, ensure_ascii=False),
                        "delta_counts": json.dumps(event.delta_counts, ensure_ascii=False),
                        "evidence_bonds": json.dumps(event.evidence_bonds, ensure_ascii=False),
                        "confidence": event.confidence,
                    }
                )

    def _write_bond_events(self, out_prefix: str, bond_events: list) -> None:
        path = f"{out_prefix}_events_bonds.csv"
        print(f"[write] {path}")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "event_id",
                    "frame",
                    "edge_type",
                    "action",
                    "i",
                    "j",
                    "distance",
                    "threshold_form",
                    "threshold_break",
                ],
            )
            writer.writeheader()
            for event in bond_events:
                writer.writerow(
                    {
                        "event_id": event.event_id,
                        "frame": event.frame,
                        "edge_type": event.edge.type.value,
                        "action": event.action,
                        "i": event.edge.i,
                        "j": event.edge.j,
                        "distance": event.evidence.distance,
                        "threshold_form": event.evidence.threshold_form,
                        "threshold_break": event.evidence.threshold_break,
                    }
                )

    def _write_main_graph(self, out_prefix: str, graph: nx.DiGraph) -> None:
        path = f"{out_prefix}_network.graphml"
        print(f"[write] {path}")
        finalize_graph_for_export(graph)
        nx.write_graphml(graph, path)

    def _build_summary_graph(self, transform_events: list) -> nx.DiGraph:
        summary_graph = ensure_graph()
        for event in transform_events:
            add_transform(
                summary_graph,
                event.reactants,
                event.products,
                event.type,
                frame=event.frame,
            )
        return summary_graph

    def _write_reaction_summary(
        self,
        out_prefix: str,
        transform_events: list,
        *,
        include_frames: bool,
    ) -> None:
        summary_graph = self._build_summary_graph(transform_events)
        lines = summarize_reversible_reactions(
            summary_graph,
            use="orig_id",
            min_total_weight=1,
            sort_by="total",
            include_frames=include_frames,
            unique_frames=False,
        )

        path = f"{out_prefix}_reactions_summary.tsv"
        print(f"[write] {path}")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))
