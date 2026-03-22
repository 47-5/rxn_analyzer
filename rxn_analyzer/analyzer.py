"""Top-level analysis orchestrator.

This module coordinates the main runtime flow:
1. build edge/bond state for the current frame
2. build the species snapshot
3. build active-site state for the frame
4. emit ordinary transform events
5. emit active-site events
6. write final outputs after the trajectory is complete

Detailed business logic should live in `species/`, `active_site/`, and `graph/`.
"""

from __future__ import annotations
from dataclasses import dataclass, field

from ase import Atoms

from .criteria import Criteria
from .tracking import EdgeTracker
from .edge_pipeline import EdgePipeline
from .species import (
    ComponentMapper,
    EventIdCounter,
    SmilesBestEffortStrategy,
    SmilesComboStrategy,
    SmilesEdgesRDKitStrategy,
    SmilesOpenBabel3DStrategy,
    SmilesRDKit3DStrategy,
    SmilesRDKitTopologyStrategy,
    SmilesStrategy,
    SpeciesLabeler,
    SpeciesPipeline,
    SpeciesRuntime,
    TransformEmitter,
)
from .frame_species_logger import FrameSpeciesLogger
from .output_writer import OutputWriter
from .baseline_gate import BaselineGate
from .sites import SiteDefinition
from .active_site import (
    ActiveSiteDefinition,
    ActiveSiteRuntime,
)

from .graph import ensure_bipartite_graph


@dataclass
class AnalyzerConfig:
    wl_iters: int = 3

    baseline_persist: int = 3
    record_init_events: bool = True
    drop_init_events: bool = False

    ads_signature_mode: str = "detailed"
    geometric_site_signature_mode: str = "none"
    geometric_site_definition: SiteDefinition | None = None
    active_site_definition: ActiveSiteDefinition | None = None

    smiles_recompute_mode: str = "always"
    smiles_mode: str = "combo"

    smiles_combo_order: list[str] = field(default_factory=lambda: ["rdkit_topology", "rdkit_3d", "openbabel_3d"])
    smiles_combo_treat_suspicious_as_failure: bool = True
    smiles_combo_prefer_charged: bool = False

    smiles_allow_charged: bool = False
    smiles_allow_dot: bool = False

    rdkit_sanitize: bool = True
    rdkit_hide_hs: bool = False
    smiles_fallback_to_formula_if_suspicious: bool = True

    smiles_strategy_impl: SmilesStrategy | None = None

    split_by_component: bool = True
    component_overlap_min: float = 0.0
    component_overlap_mode: str = "symmetric"
    require_component_mapping: bool = True
    split_component_keep_bond_evidence: bool = True

    record_frame_species: bool = False
    frame_species_include_counts: bool = True
    frame_species_include_components: bool = True
    frame_species_streaming: bool = True
    frame_species_path: str | None = None

    # NEW: whether to include forward/reverse frame lists in *_reactions_summary.tsv
    reaction_summary_include_frames: bool = True

    active_site_record_states: bool = True
    active_site_record_events: bool = True
    active_site_record_joint_reactions: bool = True
    active_site_streaming: bool = True


class ReactionAnalyzer:
    def __init__(self, criteria: Criteria, slab_def, config: AnalyzerConfig | None = None, out_prefix: str | None = None):
        self.criteria = criteria
        self.slab_def = slab_def
        self.config = config or AnalyzerConfig()
        self.out_prefix = out_prefix or "out"

        self.tracker = EdgeTracker(persist=criteria.persist, cooldown=criteria.cooldown)
        self.bond_events: list = []
        self.graph = ensure_bipartite_graph()

        self.id_counter = EventIdCounter(start=1)

        self.edge_pipeline = EdgePipeline(criteria=self.criteria, slab_def=self.slab_def, tracker=self.tracker)

        self.labeler = SpeciesLabeler(
            wl_iters=self.config.wl_iters,
            ads_signature_mode=self.config.ads_signature_mode,
            geometric_site_signature_mode=self.config.geometric_site_signature_mode,
            geometric_site_definition=self.config.geometric_site_definition,
            smiles_recompute_mode=self.config.smiles_recompute_mode,
            smiles_strategy=self._build_smiles_strategy(),
            smiles_fallback_to_formula_if_suspicious=self.config.smiles_fallback_to_formula_if_suspicious,
            smiles_allow_charged=self.config.smiles_allow_charged,
            smiles_allow_dot=self.config.smiles_allow_dot,
        )
        self.species_pipeline = SpeciesPipeline(self.labeler)

        self.mapper = ComponentMapper(
            split_by_component=bool(self.config.split_by_component),
            component_overlap_min=float(self.config.component_overlap_min),
            component_overlap_mode=str(self.config.component_overlap_mode),
            require_component_mapping=bool(self.config.require_component_mapping),
            keep_bond_evidence=bool(self.config.split_component_keep_bond_evidence),
        )

        self.emitter = TransformEmitter(self.graph, self.id_counter)
        self.species_runtime = SpeciesRuntime(mapper=self.mapper, emitter=self.emitter)

        self.baseline = BaselineGate(persist=self.config.baseline_persist)

        self.frame_logger = FrameSpeciesLogger(
            record=bool(self.config.record_frame_species),
            include_counts=bool(self.config.frame_species_include_counts),
            include_components=bool(self.config.frame_species_include_components),
            streaming=bool(self.config.frame_species_streaming),
            path=self.config.frame_species_path,
        )

        self.output_writer = OutputWriter()
        self.active_site_runtime = ActiveSiteRuntime(
            definition=self.config.active_site_definition,
            out_prefix=self.out_prefix,
            record_states=bool(self.config.active_site_record_states),
            record_events=bool(self.config.active_site_record_events),
            record_joint_reactions=bool(self.config.active_site_record_joint_reactions),
            streaming=bool(self.config.active_site_streaming),
        )

    def _build_smiles_strategy(self) -> SmilesStrategy:
        if self.config.smiles_strategy_impl is not None:
            return self.config.smiles_strategy_impl

        mode = getattr(self.config, "smiles_mode", "combo").lower().strip()

        if mode == "edges_rdkit":
            mode = "rdkit_topology"
        elif mode == "best_effort":
            mode = "combo"

        if mode == "rdkit_3d":
            return SmilesRDKit3DStrategy(hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False)))
        if mode == "openbabel_3d":
            return SmilesOpenBabel3DStrategy()
        if mode == "rdkit_topology":
            return SmilesRDKitTopologyStrategy(
                sanitize=bool(getattr(self.config, "rdkit_sanitize", True)),
                hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False)),
            )
        if mode == "combo":
            order = list(getattr(self.config, "smiles_combo_order", []))
            if not order:
                order = ["rdkit_topology", "rdkit_3d", "openbabel_3d"]

            strategies: list[SmilesStrategy] = []
            for name in order:
                nm = str(name).lower().strip()
                if nm == "rdkit_3d":
                    strategies.append(SmilesRDKit3DStrategy(hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False))))
                elif nm == "openbabel_3d":
                    strategies.append(SmilesOpenBabel3DStrategy())
                elif nm == "rdkit_topology":
                    strategies.append(
                        SmilesRDKitTopologyStrategy(
                            sanitize=bool(getattr(self.config, "rdkit_sanitize", True)),
                            hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False)),
                        )
                    )
                elif nm == "edges_rdkit":
                    strategies.append(
                        SmilesRDKitTopologyStrategy(
                            sanitize=bool(getattr(self.config, "rdkit_sanitize", True)),
                            hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False)),
                        )
                    )
                elif nm == "best_effort":
                    strategies.append(SmilesRDKit3DStrategy(hide_hs=bool(getattr(self.config, "rdkit_hide_hs", False))))
                    strategies.append(SmilesOpenBabel3DStrategy())
                else:
                    raise ValueError(f"Unknown smiles strategy name in combo order: {name}")

            return SmilesComboStrategy(
                strategies=strategies,
                treat_suspicious_as_failure=bool(getattr(self.config, "smiles_combo_treat_suspicious_as_failure", True)),
                allow_charged=bool(getattr(self.config, "smiles_allow_charged", False)),
                allow_dot=bool(getattr(self.config, "smiles_allow_dot", False)),
                prefer_charged=bool(getattr(self.config, "smiles_combo_prefer_charged", False)),
            )

        raise ValueError(f"Unknown smiles_mode: {mode}")

    def process_frame(self, frame: int, atoms: Atoms):
        edge = self.edge_pipeline.step(frame, atoms, self.id_counter)
        self.bond_events.extend(edge.bond_events)

        species_snapshot = self.species_pipeline.analyze(
            atoms,
            edge.cov_edges,
            edge.ads_edges,
            edge.slab_mask,
            edge.slab_edges,
        )
        active_batch = self.active_site_runtime.analyze_frame(frame, atoms, edge, species_snapshot)

        ads_pairs = [
            species_snapshot.component_ads_pairs.get(frozenset(component), [])
            for component in species_snapshot.components
        ]
        geometric_site_assignments = [
            species_snapshot.component_geometric_site_assignments.get(frozenset(component))
            for component in species_snapshot.components
        ]
        active_ctx = self.active_site_runtime.build_frame_context(species_snapshot.components, active_batch)
        self.frame_logger.record_frame(
            frame,
            species_snapshot.labels,
            species_snapshot.multiset,
            components=species_snapshot.components,
            ads_pairs=ads_pairs,
            geometric_site_assignments=geometric_site_assignments,
            active_site_roles=active_ctx.roles,
            active_site_ids=active_ctx.site_ids,
            active_site_summaries=active_ctx.summaries,
            is_active_site_owned=active_ctx.is_site_owned,
            is_active_site_associated=active_ctx.is_site_associated,
            active_site_memberships=active_ctx.memberships,
        )

        mode, delta = self.baseline.step(frame, species_snapshot.multiset)
        emit_reactive_init_events = (not self.config.drop_init_events) and self.config.record_init_events

        if mode == "confirmed_now":
            next_id = self.active_site_runtime.consume_frame(
                frame,
                active_batch,
                edge,
                next_event_id=self.id_counter.value,
                emit_events=emit_reactive_init_events,
            )
            self.id_counter.advance_to(next_id)
            self.species_runtime.set_previous_snapshot(species_snapshot)
            return

        if mode == "warmup":
            next_id = self.active_site_runtime.consume_frame(
                frame,
                active_batch,
                edge,
                next_event_id=self.id_counter.value,
                emit_events=emit_reactive_init_events,
            )
            self.id_counter.advance_to(next_id)
            if (delta or edge.bond_events) and (not self.config.drop_init_events) and self.config.record_init_events:
                self.species_runtime.emit_init_assignment(
                    frame=frame,
                    delta=delta,
                    bond_events=edge.bond_events,
                    add_to_graph=False,
                )
            return

        if (delta or edge.bond_events):
            self.species_runtime.emit_frame_events(
                frame=frame,
                snapshot=species_snapshot,
                delta=delta,
                bond_events=edge.bond_events,
                add_to_graph=True,
            )

        next_id = self.active_site_runtime.consume_frame(
            frame,
            active_batch,
            edge,
            next_event_id=self.id_counter.value,
            emit_events=True,
        )
        self.id_counter.advance_to(next_id)
        self.species_runtime.set_previous_snapshot(species_snapshot)

    def write_outputs(self, out_prefix: str):
        self.output_writer.write_all(
            out_prefix=out_prefix,
            transform_events=self.species_runtime.transform_events,
            bond_events=self.bond_events,
            graph=self.graph,
            frame_species_logger=self.frame_logger,
            reaction_summary_include_frames=bool(self.config.reaction_summary_include_frames),
        )
        self.active_site_runtime.write_outputs(
            out_prefix=out_prefix,
            transform_events=self.species_runtime.transform_events,
        )
