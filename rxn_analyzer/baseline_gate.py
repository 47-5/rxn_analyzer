from __future__ import annotations
from collections import Counter
from typing import Literal


def _multiset_signature(ms: Counter[str]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(ms.items()))


def _delta_multiset(prev: Counter[str], curr: Counter[str]) -> dict[str, int]:
    delta: dict[str, int] = {}
    allk = set(prev.keys()) | set(curr.keys())
    for k in allk:
        dv = curr.get(k, 0) - prev.get(k, 0)
        if dv != 0:
            delta[k] = dv
    return delta


class BaselineGate:
    """
    Detects the first stable multiset signature and gates event emission.
    """
    def __init__(self, persist: int = 3):
        self.persist = max(1, int(persist))
        self._confirmed = False
        self._candidate_sig: tuple[tuple[str, int], ...] | None = None
        self._streak = 0
        self._baseline_frame: int | None = None
        self._prev_multiset: Counter[str] | None = None

    @property
    def confirmed(self) -> bool:
        return self._confirmed

    @property
    def baseline_frame(self) -> int | None:
        return self._baseline_frame

    def step(
        self,
        frame: int,
        multiset: Counter[str],
    ) -> tuple[Literal["warmup", "confirmed_now", "normal"], dict[str, int]]:
        if self._prev_multiset is None:
            self._prev_multiset = multiset

        sig = _multiset_signature(multiset)

        if not self._confirmed:
            if self._candidate_sig is None or self._candidate_sig != sig:
                self._candidate_sig = sig
                self._streak = 1
            else:
                self._streak += 1

            if self._streak >= self.persist:
                self._confirmed = True
                self._baseline_frame = frame
                self._prev_multiset = multiset
                return "confirmed_now", {}
            else:
                delta = _delta_multiset(self._prev_multiset, multiset)
                self._prev_multiset = multiset
                return "warmup", delta

        delta = _delta_multiset(self._prev_multiset, multiset)
        self._prev_multiset = multiset
        return "normal", delta