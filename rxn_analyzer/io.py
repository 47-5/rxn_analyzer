from __future__ import annotations
from ase.io import iread


def frames(path: str, stride: int = 1):
    stride = max(1, int(stride))
    for k, atoms in enumerate(iread(path, index=":")):
        if k % stride == 0:
            yield k, atoms