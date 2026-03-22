from __future__ import annotations

from .mapping import ComponentMapper
from .emitter import TransformEmitter
from .model import SpeciesFrameSnapshot


class SpeciesRuntime:
    def __init__(self, *, mapper: ComponentMapper, emitter: TransformEmitter) -> None:
        self.mapper = mapper
        self.emitter = emitter

    @property
    def transform_events(self):
        return self.emitter.transform_events

    def set_previous_snapshot(self, snapshot: SpeciesFrameSnapshot) -> None:
        self.mapper.set_prev(snapshot.component_labels)

    def emit_init_assignment(
        self,
        *,
        frame: int,
        delta: dict[str, int],
        bond_events: list,
        add_to_graph: bool,
    ) -> None:
        self.emitter.emit_from_delta(
            frame=frame,
            delta=delta,
            bond_events=bond_events,
            ttype_override="init_assignment",
            add_to_graph=add_to_graph,
        )

    def emit_frame_events(
        self,
        *,
        frame: int,
        snapshot: SpeciesFrameSnapshot,
        delta: dict[str, int],
        bond_events: list,
        add_to_graph: bool,
    ) -> None:
        plans = self.mapper.split(snapshot.component_labels, bond_events)
        if plans:
            self.emitter.emit_from_plans(
                frame=frame,
                plans=plans,
                bond_events=bond_events,
                add_to_graph=add_to_graph,
            )
            return

        self.emitter.emit_from_delta(
            frame=frame,
            delta=delta,
            bond_events=bond_events,
            ttype_override=None,
            add_to_graph=add_to_graph,
        )
