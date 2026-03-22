from __future__ import annotations
import json
import csv
from typing import IO
from collections import Counter


class FrameSpeciesLogger:
    def __init__(
        self,
        *,
        record: bool,  # 是否启用
        include_counts: bool,  # 是否写 Counter
        include_components: bool,  # 是否写原子分量
        streaming: bool,  # 是否实时写文件
        path: str | None,  # CSV 输出路径
    ):
        self.record = record
        self.include_counts = include_counts
        self.include_components = include_components
        self.streaming = streaming
        self.path = path

        self.frame_species: list[dict] = []
        self._fh: IO[str] | None = None
        self._writer: csv.DictWriter | None = None
        self._header_written: bool = False

        if self.record and self.streaming and self.path:  # 若开启记录 + 流式写入 + 已给路径，就立即打开文件并写表头
            self._ensure_writer()

    def _ensure_writer(self, out_prefix: str | None = None) -> None:
        if not self.record or not self.streaming:  # 条件检查：未开启记录或非流式，直接返回
            return
        if self._writer is not None:  # 已存在 writer：直接返回，避免重复打开
            return

        path = self.path
        if path is None and out_prefix is not None:  # 若 self.path 为 None 且提供 out_prefix，则用 "{out_prefix}_frame_species.csv"
            path = f"{out_prefix}_frame_species.csv"
            self.path = path
        if path is None:
            return

        self._fh = open(path, "w", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(
            self._fh,
            fieldnames=[
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
            ],
        )
        self._writer.writeheader()
        self._header_written = True

        if self.frame_species:
            for r in self.frame_species:
                self._write_record(r)
            self.frame_species.clear()

        try:
            self._fh.flush()
        except Exception:
            pass

    def _write_record(self, rec: dict) -> None:
        if self._writer is None:
            return
        self._writer.writerow(  # species, counts, components 都序列化为 JSON 字符串
            {
                "frame": rec.get("frame"),
                "n_components": rec.get("n_components"),
                "species": json.dumps(rec.get("species", []), ensure_ascii=False),
                "counts": json.dumps(rec.get("counts", {}), ensure_ascii=False) if "counts" in rec else "",
                "components": json.dumps(rec.get("components", []), ensure_ascii=False) if "components" in rec else "",
                "ads_pairs": json.dumps(rec.get("ads_pairs", []), ensure_ascii=False) if "ads_pairs" in rec else "",
                "geometric_site_assignments": (
                    json.dumps(rec.get("site_assignments", []), ensure_ascii=False)
                    if "site_assignments" in rec else ""
                ),
                "active_site_roles": (
                    json.dumps(rec.get("reactive_site_roles", []), ensure_ascii=False)
                    if "reactive_site_roles" in rec else ""
                ),
                "active_site_ids": (
                    json.dumps(rec.get("reactive_site_ids", []), ensure_ascii=False)
                    if "reactive_site_ids" in rec else ""
                ),
                "active_site_summaries": (
                    json.dumps(rec.get("reactive_site_summaries", []), ensure_ascii=False)
                    if "reactive_site_summaries" in rec else ""
                ),
                "is_active_site_owned": (
                    json.dumps(rec.get("is_site_owned", []), ensure_ascii=False)
                    if "is_site_owned" in rec else ""
                ),
                "is_active_site_associated": (
                    json.dumps(rec.get("is_site_associated", []), ensure_ascii=False)
                    if "is_site_associated" in rec else ""
                ),
                "active_site_memberships": (
                    json.dumps(rec.get("reactive_site_memberships", []), ensure_ascii=False)
                    if "reactive_site_memberships" in rec else ""
                ),
            }
        )
        if self._fh is not None:
            try:
                self._fh.flush()
            except Exception:
                pass

    def record_frame(
        self,
        frame: int,
        labels: list[str],
        multiset: Counter[str],
        comps: list[list[int]],
        ads_pairs: list[list[str]] | None = None,
        site_assignments: list[dict | None] | None = None,
        reactive_site_roles: list[str] | None = None,
        reactive_site_ids: list[list[str]] | None = None,
        reactive_site_summaries: list[list[str]] | None = None,
        is_site_owned: list[bool] | None = None,
        is_site_associated: list[bool] | None = None,
        reactive_site_memberships: list[list[dict]] | None = None,
    ) -> None:
        if not self.record:
            return

        rec = {
            "frame": frame,
            "n_components": len(labels),
            "species": labels,
        }
        if self.include_counts:
            rec["counts"] = dict(multiset)
        if self.include_components:
            rec["components"] = comps
        if ads_pairs is not None:
            rec["ads_pairs"] = ads_pairs
        if site_assignments is not None:
            rec["site_assignments"] = site_assignments
        if reactive_site_roles is not None:
            rec["reactive_site_roles"] = reactive_site_roles
        if reactive_site_ids is not None:
            rec["reactive_site_ids"] = reactive_site_ids
        if reactive_site_summaries is not None:
            rec["reactive_site_summaries"] = reactive_site_summaries
        if is_site_owned is not None:
            rec["is_site_owned"] = is_site_owned
        if is_site_associated is not None:
            rec["is_site_associated"] = is_site_associated
        if reactive_site_memberships is not None:
            rec["reactive_site_memberships"] = reactive_site_memberships

        if self.streaming:
            if self._writer is None and self.path:
                self._ensure_writer()
            if self._writer is not None:
                self._write_record(rec)
                return

        self.frame_species.append(rec)

    def close(self) -> None:
        if self._fh is not None:
            try:
                self._fh.flush()
            except Exception:
                pass
            try:
                self._fh.close()
            except Exception:
                pass
        self._fh = None
        self._writer = None
        self._header_written = False

    def _write_buffer_to_file(self, path: str) -> None:  # 适用于非流式模式
        if not self.frame_species:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
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
                ],
            )
            w.writeheader()
            for r in self.frame_species:
                w.writerow(
                    {
                        "frame": r.get("frame"),
                        "n_components": r.get("n_components"),
                        "species": json.dumps(r.get("species", []), ensure_ascii=False),
                        "counts": json.dumps(r.get("counts", {}), ensure_ascii=False) if "counts" in r else "",
                        "components": json.dumps(r.get("components", []), ensure_ascii=False) if "components" in r else "",
                        "ads_pairs": json.dumps(r.get("ads_pairs", []), ensure_ascii=False) if "ads_pairs" in r else "",
                        "geometric_site_assignments": (
                            json.dumps(r.get("site_assignments", []), ensure_ascii=False)
                            if "site_assignments" in r else ""
                        ),
                        "active_site_roles": (
                            json.dumps(r.get("reactive_site_roles", []), ensure_ascii=False)
                            if "reactive_site_roles" in r else ""
                        ),
                        "active_site_ids": (
                            json.dumps(r.get("reactive_site_ids", []), ensure_ascii=False)
                            if "reactive_site_ids" in r else ""
                        ),
                        "active_site_summaries": (
                            json.dumps(r.get("reactive_site_summaries", []), ensure_ascii=False)
                            if "reactive_site_summaries" in r else ""
                        ),
                        "is_active_site_owned": (
                            json.dumps(r.get("is_site_owned", []), ensure_ascii=False)
                            if "is_site_owned" in r else ""
                        ),
                        "is_active_site_associated": (
                            json.dumps(r.get("is_site_associated", []), ensure_ascii=False)
                            if "is_site_associated" in r else ""
                        ),
                        "active_site_memberships": (
                            json.dumps(r.get("reactive_site_memberships", []), ensure_ascii=False)
                            if "reactive_site_memberships" in r else ""
                        ),
                    }
                )

    def finalize(self, out_prefix: str) -> None:
        if not self.record:
            return

        if self.streaming:
            self._ensure_writer(out_prefix=out_prefix)
            if self._writer is not None and self.frame_species:
                for r in self.frame_species:
                    self._write_record(r)
                self.frame_species.clear()
            self.close()
        else:
            path = self.path or f"{out_prefix}_frame_species.csv"
            self._write_buffer_to_file(path)
