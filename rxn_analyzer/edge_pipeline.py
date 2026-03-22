from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from ase.data import covalent_radii

from .criteria import Criteria
from .edges import Edge, EdgeEvidence, EdgeType
from .tracking import EdgeTracker, BondEvent
from .slab import SlabDefinition
from .species import EventIdCounter


def _edge(i: int, j: int, t: EdgeType) -> Edge:
    return Edge(i, j, t) if i <= j else Edge(j, i, t)


def _max_neighbor_cutoff(atoms: Atoms, criteria: Criteria) -> float:
    present = np.unique(atoms.numbers)
    rmax = float(max(covalent_radii[int(z)] for z in present))
    pair_sum_max = 2.0 * rmax
    return max(
        criteria.cov.scale_break * pair_sum_max,
        criteria.ads.scale_break * pair_sum_max,
        criteria.slab.scale_break * pair_sum_max,
    )


def _classify_edge_type(slab_mask: np.ndarray, i: int, j: int) -> EdgeType:
    si = bool(slab_mask[i])
    sj = bool(slab_mask[j])
    if si and sj:
        return EdgeType.SLAB
    if si ^ sj:
        return EdgeType.ADS
    return EdgeType.COV


def build_inst_state_strict_hysteresis(
    atoms: Atoms,
    slab_mask: np.ndarray,
    criteria: Criteria,
    tracker: EdgeTracker,
) -> tuple[dict[Edge, tuple[int, EdgeEvidence]], set[Edge]]:
    """
    Strict hysteresis instant state:
      if confirmed OFF: ON when d <= d_form else OFF
      if confirmed ON : ON when d <= d_break else OFF

    Also applies optional max_neighbors pruning on *instantaneous ON edges* per type.
    """
    numbers = atoms.numbers
    max_cut = _max_neighbor_cutoff(atoms, criteria)
    i_list, j_list, d_list = neighbor_list("ijd", atoms, cutoff=max_cut)

    inst: dict[Edge, tuple[int, EdgeEvidence]] = {}
    universe: set[Edge] = set()

    edges_on_by_type: dict[EdgeType, list[tuple[Edge, float]]] = {t: [] for t in EdgeType}

    for i, j, d in zip(i_list, j_list, d_list):
        if i == j:
            continue
        et = _classify_edge_type(slab_mask, int(i), int(j))
        params = criteria.cov if et == EdgeType.COV else criteria.ads if et == EdgeType.ADS else criteria.slab

        e = _edge(int(i), int(j), et)
        universe.add(e)

        Zi = int(numbers[e.i])
        Zj = int(numbers[e.j])
        d_form, d_break = params.thresholds(Zi, Zj)

        confirmed = tracker.get_confirmed(e)

        if confirmed == 0:
            on = 1 if float(d) <= d_form else 0
        else:
            on = 1 if float(d) <= d_break else 0

        ev = EdgeEvidence(distance=float(d), threshold_form=float(d_form), threshold_break=float(d_break))
        inst[e] = (on, ev)
        if on:
            edges_on_by_type[et].append((e, float(d)))

    def prune(et: EdgeType, max_deg: int | None):
        if max_deg is None:
            return
        inc: dict[int, list[tuple[float, Edge]]] = {}
        for e, dist in edges_on_by_type[et]:
            inc.setdefault(e.i, []).append((dist, e))
            inc.setdefault(e.j, []).append((dist, e))
        keep = set()
        for _u, lst in inc.items():
            for dist, e in sorted(lst, key=lambda x: x[0])[:max_deg]:
                keep.add(e)
        for e, _dist in edges_on_by_type[et]:
            if e not in keep:
                old_on, ev = inst[e]
                if old_on == 1:
                    inst[e] = (0, ev)

    prune(EdgeType.SLAB, criteria.slab.max_neighbors)
    prune(EdgeType.ADS, criteria.ads.max_neighbors)
    prune(EdgeType.COV, criteria.cov.max_neighbors)

    return inst, universe


@dataclass
class EdgeStep:
    slab_mask: np.ndarray
    slab_edges: set[tuple[int, int]]
    cov_edges: set[tuple[int, int]]
    ads_edges: set[tuple[int, int]]
    bond_events: list[BondEvent]


class EdgePipeline:
    def __init__(self, criteria: Criteria, slab_def: SlabDefinition, tracker: EdgeTracker):
        self.criteria = criteria
        self.slab_def = slab_def
        self.tracker = tracker

    def step(self, frame: int, atoms: Atoms, id_counter: EventIdCounter) -> EdgeStep:
        slab_mask = self.slab_def.mask(atoms)
        inst, universe = build_inst_state_strict_hysteresis(atoms, slab_mask, self.criteria, self.tracker)

        bond_events_frame, next_id = self.tracker.update(
            frame=frame,
            inst_state=inst,
            universe=universe,
            next_event_id=id_counter.value,
        )
        id_counter.advance_to(next_id)

        cov_edges = self.tracker.confirmed_edges(EdgeType.COV)
        ads_edges = self.tracker.confirmed_edges(EdgeType.ADS)
        slab_edges = self.tracker.confirmed_edges(EdgeType.SLAB)

        return EdgeStep(
            slab_mask=slab_mask,
            slab_edges=slab_edges,
            cov_edges=cov_edges,
            ads_edges=ads_edges,
            bond_events=bond_events_frame,
        )
