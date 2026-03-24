from __future__ import annotations

import csv
import json
from collections import defaultdict

import networkx as nx

from .graph import finalize_graph_for_export


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

    def _summarize_reversible_reactions_from_events(
        self,
        transform_events: list,
        *,
        include_frames: bool,
    ) -> list[str]:
        pair_counts: dict[tuple[str, str], int] = defaultdict(int)
        pair_frames: dict[tuple[str, str], list[int]] = defaultdict(list)

        for event in transform_events:
            reactants = " + ".join(sorted(event.reactants))
            products = " + ".join(sorted(event.products))
            pair_counts[(reactants, products)] += 1
            pair_frames[(reactants, products)].append(int(event.frame))

        seen: set[frozenset[str]] = set()
        rows: list[dict[str, object]] = []
        for reactants, products in pair_counts:
            pair = frozenset((reactants, products))
            if pair in seen:
                continue
            seen.add(pair)

            if (reactants, products) <= (products, reactants):
                lhs, rhs = reactants, products
            else:
                lhs, rhs = products, reactants

            fwd = pair_counts.get((lhs, rhs), 0)
            rev = pair_counts.get((rhs, lhs), 0)
            total = fwd + rev
            rows.append(
                {
                    "lhs": lhs,
                    "rhs": rhs,
                    "fwd": fwd,
                    "rev": rev,
                    "total": total,
                    "fwd_frames": pair_frames.get((lhs, rhs), []),
                    "rev_frames": pair_frames.get((rhs, lhs), []),
                }
            )

        rows.sort(key=lambda row: int(row["total"]), reverse=True)
        if include_frames:
            lines = ["reaction\tfwd\trev\ttotal\tfwd_frames\trev_frames"]
            for row in rows:
                lines.append(
                    f"{row['lhs']}<->{row['rhs']}\t{row['fwd']}\t{row['rev']}\t{row['total']}\t"
                    f"{json.dumps(row['fwd_frames'], ensure_ascii=False)}\t"
                    f"{json.dumps(row['rev_frames'], ensure_ascii=False)}"
                )
            return lines

        lines = ["reaction\tfwd\trev\ttotal"]
        for row in rows:
            lines.append(f"{row['lhs']}<->{row['rhs']}\t{row['fwd']}\t{row['rev']}\t{row['total']}")
        return lines

    def _write_reaction_summary(
        self,
        out_prefix: str,
        transform_events: list,
        *,
        include_frames: bool,
    ) -> None:
        lines = self._summarize_reversible_reactions_from_events(
            transform_events,
            include_frames=include_frames,
        )

        path = f"{out_prefix}_reactions_summary.tsv"
        print(f"[write] {path}")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))
