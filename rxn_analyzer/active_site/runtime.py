"""Runtime facade for dynamic active-site analysis.

This layer connects:
- single-frame active-site state recognition
- active-site membership summaries for frame-level logging
- cross-frame active-site event tracking
- optional streaming/final active-site outputs
"""

from __future__ import annotations

from ase import Atoms

from .integration import build_active_site_memberships, summarize_active_site_memberships
from .model import ActiveSiteDefinition, ActiveSiteStateFrame
from .output import ActiveSiteOutputWriter
from .pipeline import ActiveSitePipeline
from .tracker import ActiveSiteTracker
from .types import ActiveSiteFrameBatch, ActiveSiteFrameContext
from ..species import SpeciesFrameSnapshot


class ActiveSiteRuntime:
    def __init__(
        self,
        *,
        definition: ActiveSiteDefinition | None,
        out_prefix: str,
        record_states: bool,
        record_events: bool,
        record_joint_reactions: bool,
        streaming: bool,
    ) -> None:
        self.definition = definition
        self.out_prefix = out_prefix
        self.record_states = bool(record_states)
        self.record_events = bool(record_events)
        self.record_joint_reactions = bool(record_joint_reactions)
        self.streaming = bool(streaming)

        self.pipeline = ActiveSitePipeline(definition) if definition is not None else None
        self.tracker = ActiveSiteTracker()
        self.output_writer = ActiveSiteOutputWriter()
        self.state_frames: list[ActiveSiteStateFrame] = []

    @property
    def enabled(self) -> bool:
        return self.pipeline is not None

    def analyze_frame(
        self,
        frame: int,
        atoms: Atoms,
        edge,
        species_snapshot: SpeciesFrameSnapshot,
    ) -> ActiveSiteFrameBatch | None:
        if self.pipeline is None:
            return None
        return self.pipeline.analyze(
            frame=frame,
            atoms=atoms,
            slab_mask=edge.slab_mask,
            cov_edges=edge.cov_edges,
            ads_edges=edge.ads_edges,
            components=species_snapshot.components,
            labels=species_snapshot.labels,
        )

    def build_frame_context(
        self,
        comps: list[list[int]],
        batch: ActiveSiteFrameBatch | None,
    ) -> ActiveSiteFrameContext:
        states = [] if batch is None else batch.states
        memberships = build_active_site_memberships(comps, states)
        roles, site_ids, summaries, is_site_owned, is_site_associated = summarize_active_site_memberships(memberships)
        return ActiveSiteFrameContext(
            memberships=memberships,
            roles=roles,
            site_ids=site_ids,
            summaries=summaries,
            is_site_owned=is_site_owned,
            is_site_associated=is_site_associated,
        )

    def consume_frame(self, frame: int, batch: ActiveSiteFrameBatch | None, edge, *, next_event_id: int, emit_events: bool) -> int:
        if batch is None:
            return next_event_id

        if self.record_states:
            self.state_frames.extend(batch.states)
            if self.streaming:
                self.output_writer.record_state_frames(self.out_prefix, batch.states)

        prev_event_count = len(self.tracker.events)
        prev_joint_count = len(self.tracker.joint_reactions)
        next_event_id = self.tracker.update(
            frame=frame,
            states=batch.states,
            bond_events=edge.bond_events,
            next_event_id=next_event_id,
            emit_events=emit_events,
        )

        if self.streaming:
            if self.record_events:
                self.output_writer.record_events(self.out_prefix, self.tracker.events[prev_event_count:])
            if self.record_joint_reactions:
                self.output_writer.record_joint_reactions(
                    self.out_prefix,
                    self.tracker.joint_reactions[prev_joint_count:],
                )
        return next_event_id

    def write_outputs(self, *, out_prefix: str, transform_events: list) -> None:
        if not self.enabled:
            return
        self.output_writer.write_all(
            out_prefix=out_prefix,
            transform_events=transform_events,
            state_frames=self.state_frames if self.record_states else [],
            events=self.tracker.events if self.record_events else [],
            joint_reactions=self.tracker.joint_reactions if self.record_joint_reactions else [],
            streaming=self.streaming,
        )
