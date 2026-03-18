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
            fieldnames=["frame", "n_components", "species", "counts", "components"],
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
                fieldnames=["frame", "n_components", "species", "counts", "components"],
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