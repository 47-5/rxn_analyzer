# Architecture Notes

`rxn_analyzer` currently has four core runtime layers:

1. `config_loader`
   Reads YAML, validates new-style config keys, builds runtime objects.
2. `analyzer`
   Orchestrates the full per-frame workflow.
3. `species`
   Handles non-host connected components, labeling, cross-frame mapping, and transform emission.
4. `active_site`
   Handles dynamic active-site states, events, and coupling to ordinary reactions.
5. `graph`
   Handles reaction-graph node IDs, graph attributes, and GraphML-ready reaction-network construction.

## Main Data Flow

Per frame, the runtime is roughly:

1. Trajectory frame enters [`ReactionAnalyzer`](/D:/code/rxn_analyzer/rxn_analyzer/analyzer.py).
2. [`EdgePipeline`](/D:/code/rxn_analyzer/rxn_analyzer/edge_pipeline.py) builds:
   - covalent edges
   - adsorption edges
   - host/host edges
   - bond events
3. [`SpeciesPipeline`](/D:/code/rxn_analyzer/rxn_analyzer/species/pipeline.py) builds a [`SpeciesFrameSnapshot`](/D:/code/rxn_analyzer/rxn_analyzer/species/model.py).
4. [`SpeciesRuntime`](/D:/code/rxn_analyzer/rxn_analyzer/species/runtime.py) compares current and previous snapshots and emits transform events.
5. [`ActiveSiteRuntime`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/runtime.py) computes current active-site states and cross-frame active-site events.
6. [`FrameSpeciesLogger`](/D:/code/rxn_analyzer/rxn_analyzer/frame_species_logger.py) records optional per-frame CSV output.
7. [`OutputWriter`](/D:/code/rxn_analyzer/rxn_analyzer/output_writer.py) and [`ActiveSiteOutputWriter`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/output.py) write final CSV/GraphML outputs.

## Module Boundaries

### `rxn_analyzer/analyzer.py`

Keep this file as the orchestration layer.
It should coordinate pipelines and runtimes, but avoid owning detailed business logic.

### `rxn_analyzer/species/`

- `chem.py`: single-frame chemistry helpers
- `model.py`: `SpeciesFrameSnapshot`
- `pipeline.py`: single-frame species identification
- `mapping.py`: cross-frame component matching
- `emitter.py`: transform-event emission
- `runtime.py`: species runtime facade used by `analyzer.py`

### `rxn_analyzer/active_site/`

- `model.py`: public active-site data models
- `types.py`: intermediate runtime types
- `pipeline.py`: single-frame state analysis
- `tracker.py`: cross-frame active-site event tracking
- `coupling.py`: couple active-site state changes with ordinary transform events
- `runtime.py`: active-site runtime facade used by `analyzer.py`
- `output.py`: active-site CSV outputs and site-aware graph export

### `rxn_analyzer/graph/`

- `ids.py`: graph node/reaction ID allocation
- `attrs.py`: graph node/edge attribute helpers
- `network.py`: reaction-graph construction and reversible-summary logic

## Naming Conventions

The current codebase intentionally uses:

- `host_definition`
- `geometric_site`
- `active_site`
- `AnalyzerConfig.geometric_site_*`
- `AnalyzerConfig.active_site_*`

Legacy config keys such as `slab_definition`, `site`, and `reactive_site` are no longer accepted.

## Practical Reading Order

If you are debugging the main analysis path, read in this order:

1. [`rxn_analyzer/analyzer.py`](/D:/code/rxn_analyzer/rxn_analyzer/analyzer.py)
2. [`rxn_analyzer/species/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/runtime.py)
3. [`rxn_analyzer/active_site/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/runtime.py)
4. [`rxn_analyzer/graph/network.py`](/D:/code/rxn_analyzer/rxn_analyzer/graph/network.py)
5. Then drill into `pipeline / mapping / emitter / tracker / coupling` as needed.
