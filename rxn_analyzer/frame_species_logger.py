from __future__ import annotations

import csv
import json
from collections import Counter
from typing import IO


class FrameSpeciesLogger:
    FIELDNAMES = [
        "frame",
        "n_components",
        "species",
        "counts",
        "components",
        "ads_pairs",
        "geometric_site_assignments",
        "active_site_roles",
        "active_site_ids",
        "active_site_summaries",
        "is_active_site_owned",
        "is_active_site_associated",
        "active_site_memberships",
    ]

    def __init__(
        self,
        *,
        record: bool,
        include_counts: bool,
        include_components: bool,
        streaming: bool,
        path: str | None,
    ) -> None:
        self.record = record
        self.include_counts = include_counts
        self.include_components = include_components
        self.streaming = streaming
        self.path = path

        self.frame_species: list[dict] = []
        self._fh: IO[str] | None = None
        self._writer: csv.DictWriter | None = None

        if self.record and self.streaming and self.path:
            self._ensure_writer()

    def _ensure_writer(self, out_prefix: str | None = None) -> None:
        if not self.record or not self.streaming:
            return
        if self._writer is not None:
            return

        path = self.path
        if path is None and out_prefix is not None:
            path = f"{out_prefix}_frame_species.csv"
            self.path = path
        if path is None:
            return

        self._fh = open(path, "w", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(self._fh, fieldnames=self.FIELDNAMES)
        self._writer.writeheader()

        if self.frame_species:
            for record in self.frame_species:
                self._write_record(record)
            self.frame_species.clear()

        self._flush()

    def _flush(self) -> None:
        if self._fh is None:
            return
        try:
            self._fh.flush()
        except Exception:
            pass

    @staticmethod
    def _json_or_empty(record: dict, key: str, default) -> str:
        return json.dumps(record.get(key, default), ensure_ascii=False) if key in record else ""

    def _csv_row_from_record(self, record: dict) -> dict:
        return {
            "frame": record.get("frame"),
            "n_components": record.get("n_components"),
            "species": json.dumps(record.get("species", []), ensure_ascii=False),
            "counts": self._json_or_empty(record, "counts", {}),
            "components": self._json_or_empty(record, "components", []),
            "ads_pairs": self._json_or_empty(record, "ads_pairs", []),
            "geometric_site_assignments": self._json_or_empty(record, "geometric_site_assignments", []),
            "active_site_roles": self._json_or_empty(record, "active_site_roles", []),
            "active_site_ids": self._json_or_empty(record, "active_site_ids", []),
            "active_site_summaries": self._json_or_empty(record, "active_site_summaries", []),
            "is_active_site_owned": self._json_or_empty(record, "is_active_site_owned", []),
            "is_active_site_associated": self._json_or_empty(record, "is_active_site_associated", []),
            "active_site_memberships": self._json_or_empty(record, "active_site_memberships", []),
        }

    def _write_record(self, record: dict) -> None:
        if self._writer is None:
            return
        self._writer.writerow(self._csv_row_from_record(record))
        self._flush()

    def record_frame(
        self,
        frame: int,
        labels: list[str],
        multiset: Counter[str],
        components: list[list[int]],
        ads_pairs: list[list[str]] | None = None,
        geometric_site_assignments: list[dict | None] | None = None,
        active_site_roles: list[str] | None = None,
        active_site_ids: list[list[str]] | None = None,
        active_site_summaries: list[list[str]] | None = None,
        is_active_site_owned: list[bool] | None = None,
        is_active_site_associated: list[bool] | None = None,
        active_site_memberships: list[list[dict]] | None = None,
    ) -> None:
        if not self.record:
            return

        record = {
            "frame": frame,
            "n_components": len(labels),
            "species": labels,
        }
        if self.include_counts:
            record["counts"] = dict(multiset)
        if self.include_components:
            record["components"] = components
        if ads_pairs is not None:
            record["ads_pairs"] = ads_pairs
        if geometric_site_assignments is not None:
            record["geometric_site_assignments"] = geometric_site_assignments
        if active_site_roles is not None:
            record["active_site_roles"] = active_site_roles
        if active_site_ids is not None:
            record["active_site_ids"] = active_site_ids
        if active_site_summaries is not None:
            record["active_site_summaries"] = active_site_summaries
        if is_active_site_owned is not None:
            record["is_active_site_owned"] = is_active_site_owned
        if is_active_site_associated is not None:
            record["is_active_site_associated"] = is_active_site_associated
        if active_site_memberships is not None:
            record["active_site_memberships"] = active_site_memberships

        if self.streaming:
            if self._writer is None and self.path:
                self._ensure_writer()
            if self._writer is not None:
                self._write_record(record)
                return

        self.frame_species.append(record)

    def close(self) -> None:
        if self._fh is not None:
            self._flush()
            try:
                self._fh.close()
            except Exception:
                pass
        self._fh = None
        self._writer = None

    def _write_buffer_to_file(self, path: str) -> None:
        if not self.frame_species:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=self.FIELDNAMES)
            writer.writeheader()
            for record in self.frame_species:
                writer.writerow(self._csv_row_from_record(record))

    def finalize(self, out_prefix: str) -> None:
        if not self.record:
            return

        if self.streaming:
            self._ensure_writer(out_prefix=out_prefix)
            if self._writer is not None and self.frame_species:
                for record in self.frame_species:
                    self._write_record(record)
                self.frame_species.clear()
            self.close()
            return

        path = self.path or f"{out_prefix}_frame_species.csv"
        self._write_buffer_to_file(path)
