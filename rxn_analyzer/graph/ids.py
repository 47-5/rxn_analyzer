from __future__ import annotations

_NODE_ID_MAP: dict[str, str] = {}
_NEXT_NODE_ID: int = 0
_RXN_NODE_ID_MAP: dict[tuple, str] = {}
_NEXT_RXN_ID: int = 0


def reset_node_id_mapping() -> None:
    global _NODE_ID_MAP, _NEXT_NODE_ID, _RXN_NODE_ID_MAP, _NEXT_RXN_ID
    _NODE_ID_MAP = {}
    _NEXT_NODE_ID = 0
    _RXN_NODE_ID_MAP = {}
    _NEXT_RXN_ID = 0


def get_or_create_numeric_node_id(long_id: str) -> str:
    global _NEXT_NODE_ID
    if long_id in _NODE_ID_MAP:
        return _NODE_ID_MAP[long_id]
    num_id = str(_NEXT_NODE_ID)
    _NODE_ID_MAP[long_id] = num_id
    _NEXT_NODE_ID += 1
    return num_id


def reaction_signature(reactants_sorted: list[str], products_sorted: list[str], ttype: str) -> tuple:
    return (tuple(reactants_sorted), tuple(products_sorted), str(ttype or ""))


def get_or_create_reaction_node_id(sig: tuple) -> str:
    global _NEXT_RXN_ID
    if sig in _RXN_NODE_ID_MAP:
        return _RXN_NODE_ID_MAP[sig]
    rid = f"r{_NEXT_RXN_ID}"
    _RXN_NODE_ID_MAP[sig] = rid
    _NEXT_RXN_ID += 1
    return rid
