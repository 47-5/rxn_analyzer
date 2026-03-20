from __future__ import annotations
from dataclasses import dataclass, field

from ase import Atoms

from .criteria import Criteria
from .tracking import EdgeTracker
from .edge_pipeline import EdgePipeline
from .species_pipeline import SpeciesLabeler, SpeciesPipeline
from .component_mapper import ComponentMapper
from .event_emitter import EventIdCounter, TransformEmitter
from .frame_species_logger import FrameSpeciesLogger
from .output_writer import OutputWriter
from .baseline_gate import BaselineGate
from .sites import SiteDefinition
from .site_model import ReactiveSiteDefinition, ReactiveSiteStateFrame
from .site_pipeline import ReactiveSitePipeline
from .site_tracker import ReactiveSiteTracker
from .site_output import ReactiveSiteOutputWriter

from .species import (
    SmilesStrategy,
    SmilesBestEffortStrategy,
    SmilesEdgesRDKitStrategy,
    SmilesRDKit3DStrategy,
    SmilesOpenBabel3DStrategy,
    SmilesRDKitTopologyStrategy,
    SmilesComboStrategy,
)
from .network import ensure_bipartite_graph
from .events import TransformEvent


@dataclass
class AnalyzerConfig:
    wl_iters: int = 3

    baseline_persist: int = 3
    record_init_events: bool = True
    drop_init_events: bool = False

    ads_signature_mode: str = "detailed"
    site_signature_mode: str = "none"
    site_definition: SiteDefinition | None = None
    reactive_site_definition: ReactiveSiteDefinition | None = None

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

    reactive_site_record_states: bool = True
    reactive_site_record_events: bool = True
    reactive_site_record_joint_reactions: bool = True


class ReactionAnalyzer:
    def __init__(self, criteria: Criteria, slab_def, config: AnalyzerConfig | None = None):
        self.criteria = criteria
        self.slab_def = slab_def
        self.config = config or AnalyzerConfig()

        self.tracker = EdgeTracker(persist=criteria.persist, cooldown=criteria.cooldown)
        self.bond_events: list = []
        self.graph = ensure_bipartite_graph()

        self.id_counter = EventIdCounter(start=1)

        self.edge_pipeline = EdgePipeline(criteria=self.criteria, slab_def=self.slab_def, tracker=self.tracker)

        self.labeler = SpeciesLabeler(
            wl_iters=self.config.wl_iters,
            ads_signature_mode=self.config.ads_signature_mode,
            site_signature_mode=self.config.site_signature_mode,
            site_definition=self.config.site_definition,
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

        self.baseline = BaselineGate(persist=self.config.baseline_persist)

        self.frame_logger = FrameSpeciesLogger(
            record=bool(self.config.record_frame_species),
            include_counts=bool(self.config.frame_species_include_counts),
            include_components=bool(self.config.frame_species_include_components),
            streaming=bool(self.config.frame_species_streaming),
            path=self.config.frame_species_path,
        )

        self.output_writer = OutputWriter()
        self.reactive_site_pipeline = (
            ReactiveSitePipeline(self.config.reactive_site_definition)
            if self.config.reactive_site_definition is not None
            else None
        )
        self.reactive_site_tracker = ReactiveSiteTracker()
        self.reactive_site_output_writer = ReactiveSiteOutputWriter()
        self.reactive_site_states: list[ReactiveSiteStateFrame] = []

    def _process_reactive_sites(self, frame: int, atoms: Atoms, edge, sp) -> None:
        if self.reactive_site_pipeline is None:
            return

        batch = self.reactive_site_pipeline.analyze(
            frame=frame,
            atoms=atoms,
            slab_mask=edge.slab_mask,
            cov_edges=edge.cov_edges,
            ads_edges=edge.ads_edges,
            components=sp.comps,
            labels=sp.labels,
        )
        if self.config.reactive_site_record_states:
            self.reactive_site_states.extend(batch.states)
        next_id = self.reactive_site_tracker.update(
            frame=frame,
            states=batch.states,
            bond_events=edge.bond_events,
            next_event_id=self.id_counter.value,
        )
        self.id_counter.advance_to(next_id)

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

        sp = self.species_pipeline.analyze(atoms, edge.cov_edges, edge.ads_edges, edge.slab_mask)

        self.frame_logger.record_frame(frame, sp.labels, sp.multiset, sp.comps)

        mode, delta = self.baseline.step(frame, sp.multiset)

        if mode == "confirmed_now":
            self._process_reactive_sites(frame, atoms, edge, sp)
            self.mapper.set_prev(sp.comp_labels)
            return

        if mode == "warmup":
            self._process_reactive_sites(frame, atoms, edge, sp)
            if (delta or edge.bond_events) and (not self.config.drop_init_events) and self.config.record_init_events:
                self.emitter.emit_from_delta(
                    frame=frame,
                    delta=delta,
                    bond_events=edge.bond_events,
                    ttype_override="init_assignment",
                    add_to_graph=False,
                )
            return

        if (delta or edge.bond_events):
            plans = self.mapper.split(sp.comp_labels, edge.bond_events)
            if plans:
                self.emitter.emit_from_plans(frame, plans, edge.bond_events, add_to_graph=True)
            else:
                self.emitter.emit_from_delta(
                    frame=frame,
                    delta=delta,
                    bond_events=edge.bond_events,
                    ttype_override=None,
                    add_to_graph=True,
                )

        self._process_reactive_sites(frame, atoms, edge, sp)
        self.mapper.set_prev(sp.comp_labels)

    def write_outputs(self, out_prefix: str):
        self.output_writer.write_all(
            out_prefix=out_prefix,
            transform_events=self.emitter.transform_events,
            bond_events=self.bond_events,
            graph=self.graph,
            frame_species_logger=self.frame_logger,
            reaction_summary_include_frames=bool(self.config.reaction_summary_include_frames),
        )
        if self.reactive_site_pipeline is not None:
            self.reactive_site_output_writer.write_all(
                out_prefix=out_prefix,
                transform_events=self.emitter.transform_events,
                state_frames=self.reactive_site_states if self.config.reactive_site_record_states else [],
                events=self.reactive_site_tracker.events if self.config.reactive_site_record_events else [],
                joint_reactions=(
                    self.reactive_site_tracker.joint_reactions
                    if self.config.reactive_site_record_joint_reactions
                    else []
                ),
            )
