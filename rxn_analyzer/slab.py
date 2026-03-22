from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class SlabDefinition:
    """
    Host/slab atom selection.

    Internal atom indices are always normalized to 0-based before they are stored.
    If both indices and elements are provided, indices take precedence.
    """

    indices: frozenset[int] | None = None
    elements: frozenset[str] | None = None
    invert: bool = False

    def mask(self, atoms: Atoms) -> np.ndarray:
        n = len(atoms)
        m = np.zeros(n, dtype=bool)

        if self.indices is not None:
            idx = [i for i in self.indices if 0 <= i < n]
            m[idx] = True
            return ~m if self.invert else m

        if self.elements is not None:
            s = set(self.elements)
            syms = atoms.get_chemical_symbols()
            for i, e in enumerate(syms):
                if e in s:
                    m[i] = True

        return ~m if self.invert else m


# Preferred semantic alias: host/background framework definition.
HostDefinition = SlabDefinition
