from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class SlabDefinition:
    """
    只要 indices 不为 None，elements 就被忽略
    """
    indices: frozenset[int] | None = None  # 定义 slab 原子的识别方式，有两种互斥方式
    elements: frozenset[str] | None = None  # 定义 slab 原子的识别方式，有两种互斥方式

    def mask(self, atoms: Atoms) -> np.ndarray:
        n = len(atoms)
        m = np.zeros(n, dtype=bool)

        if self.indices is not None:
            idx = [i for i in self.indices if 0 <= i < n]
            m[idx] = True
            return m

        if self.elements is not None:
            s = set(self.elements)
            syms = atoms.get_chemical_symbols()
            for i, e in enumerate(syms):
                if e in s:
                    m[i] = True

        return m