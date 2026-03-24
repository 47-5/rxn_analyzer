# rxn_analyzer User Manual

## 1. Purpose of This Manual

This manual explains the design goals, architecture, core concepts, usage patterns, outputs, suitable application scenarios, and graph postprocessing workflow of `rxn_analyzer`.

The program is designed for researchers who want to turn atomistic trajectories into chemically interpretable events and reaction networks. It is especially useful when the main challenge is not reading coordinates themselves, but answering questions such as:

- Which bonds were formed or broken?
- When did a specific intermediate appear?
- Is a species in the gas phase, adsorbed, or part of a site state?
- Did a reaction occur on a specific site?
- Did the active site itself change during the reaction?

## 2. Design Goals

The core design idea of `rxn_analyzer` is to convert continuous trajectory data into structured chemical information in three steps:

1. Recover chemically meaningful connectivity from geometry.
2. Recognize framewise species and site states.
3. Compare neighboring frames to extract events, transitions, and graphs.

The software is not just a generic MD postprocessor. It is designed around:

- reaction-event extraction
- reaction-network construction
- geometric-site annotation
- active-site state tracking and coupling analysis

## 3. Core Concepts

### 3.1 host

`host` represents the structural background of the system, such as:

- a metal slab
- a zeolite framework
- a support or extended solid host

Its main role is to define what belongs to the host and what belongs to non-host chemistry.

### 3.2 species

`species` are connected non-host components recognized frame by frame. These may correspond to:

- gas-phase molecules
- adsorbed intermediates
- products
- any non-host connected molecular fragment

### 3.3 geometric_site

`geometric_site` is a static geometric annotation layer. It is intended for cases such as:

- top / bridge / hollow sites on metal surfaces
- fixed pore-mouth templates
- other predefined geometric adsorption positions

Its main purpose is to answer: where is a species adsorbed?

### 3.4 active_site

`active_site` is a dynamic chemically meaningful center. It is intended for cases such as:

- BAS / LAS in zeolite systems
- local metal centers such as `Ga(OH)2`
- defect-associated reactive centers

Each active site contains:

- `core_members`: stable host atoms defining site identity
- `state_members`: dynamic atoms defining current site state

Within `state_members`, the software distinguishes:

- `intrinsic`
- `incorporated`
- `associated`

This allows the program to distinguish between ?something is adsorbed on the site? and ?something has become part of the site state?.

## 4. High-Level Architecture

The main runtime layers are:

1. `config_loader`
2. `analyzer`
3. `species`
4. `active_site`
5. `graph`

### 4.1 config_loader

`config_loader` reads YAML files, validates the supported keys, and builds runtime configuration objects.

The current configuration naming is standardized to:

- `host_definition`
- `geometric_site`
- `active_site`

### 4.2 analyzer

`analyzer` is the orchestration layer. It coordinates:

- edge recognition
- species analysis
- active-site analysis
- event extraction
- output writing

### 4.3 species submodule

The `species/` submodule handles the normal non-host reaction line. Its main files are:

- `chem.py`
- `model.py`
- `pipeline.py`
- `mapping.py`
- `emitter.py`
- `runtime.py`

### 4.4 active_site submodule

The `active_site/` submodule handles active-center chemistry. Its main files are:

- `model.py`
- `types.py`
- `pipeline.py`
- `tracker.py`
- `coupling.py`
- `runtime.py`
- `output.py`
- `rules.py`

### 4.5 graph submodule

The `graph/` submodule handles GraphML construction and graph postprocessing.

## 5. Framewise Analysis Workflow

Each frame is processed through the following stages:

1. `EdgePipeline` identifies covalent, adsorption, and host-host edges.
2. The species pipeline builds a `SpeciesFrameSnapshot`.
3. The species runtime compares current and previous snapshots to emit transform events.
4. The active-site runtime recognizes site states and site-related events.
5. The output layer writes CSV files and GraphML graphs.

## 6. Configuration Overview

The main YAML sections are:

- `run`
- `criteria`
- `host_definition`
- `geometric_site`
- `active_site`
- `analyzer`

### 6.1 host_definition

`host_definition` supports:

- element-based selection
- index-based selection
- inverted selection
- index ranges such as `"1-20"`
- `0-based` or `1-based` indexing

### 6.2 geometric_site

`geometric_site` is defined in a separate YAML file and can optionally be written into species labels through `analyzer.geometric_site_signature_mode`.

### 6.3 active_site

`active_site` is also defined in a separate YAML file. A typical BAS-like definition may look like:

```yaml
- id: "B_Ga_01"
  family: "BAS_Ga"
  core_members: [616, 31, 35]
  initial_state_members: [608]
  allowed_state_elements: ["H", "O", "C"]
  max_state_members: 20
  rule_profile: "bas_default"
```

This pattern is especially useful in zeolite and local-acid-center systems.

## 7. Main Outputs

The main outputs include:

- `<out_prefix>_events_transform.csv`
- `<out_prefix>_events_bonds.csv`
- `<out_prefix>_network.graphml`
- `<out_prefix>_reactions_summary.tsv`
- `<out_prefix>_frame_species.csv`

With `active_site` enabled, additional outputs include:

- `<out_prefix>_site_states.csv`
- `<out_prefix>_site_events.csv`
- `<out_prefix>_joint_site_reactions.csv`
- `<out_prefix>_site_reaction_coupling.csv`
- `<out_prefix>_active_site_state_counts.csv`
- `<out_prefix>_active_site_family_state_counts.csv`
- `<out_prefix>_active_site_transitions.csv`
- `<out_prefix>_active_site_family_transitions.csv`
- `<out_prefix>_reaction_by_active_site_state.csv`
- `<out_prefix>_reaction_by_active_site_family_state.csv`
- `<out_prefix>_network_site_aware.graphml`
- `<out_prefix>_site_outputs_README.txt`

## 8. Active-Site State Semantics

The current active-site rules already support more chemical state labels than simple formulas. For example, metal-centered LAS-like sites may now produce labels such as:

- `ga_h`
- `ga_h2`
- `ga_oh`
- `ga_oh2`
- `ga_h_oh`

BAS-like sites can produce labels such as:

- `protonated_bas`
- `deprotonated_bas`
- `alkoxy_bas`

This makes the site-state tables and site-aware graph easier to interpret chemically.

## 9. Site-Aware Interpretation

The software is designed to help distinguish three important cases:

- reaction occurs while site state does not change
- reaction occurs and site state changes
- site state changes without a clear normal reaction event

This distinction is central when studying whether a site acts only as an environment or participates directly in the chemistry.

## 10. Graph Postprocessing

Graph postprocessing is now task-oriented and YAML-driven. The three supported modes are:

- `focus`
- `context`
- `story`

### 10.1 focus

`focus` extracts the most important backbone from a full bipartite graph. Key options include:

- `score_mode`
- `top_species`
- `top_reactions`
- `top_species_per_family`
- `species_family_mode`
- `min_edge_weight`
- `collapse_reversible`
- `include_neighbors`
- `keep_nodes`
- `protect_seed_neighbors`

The current `species_family_mode` implementation supports `ads_state`, which splits species into `ads` and `non_ads`.

### 10.2 context

`context` extracts a local subgraph around one or more seeds. Seeds can be matched by:

- node id
- `orig_id`
- `label`
- `display_label`

### 10.3 story

`story` extracts source-to-target path subgraphs. It supports:

- `path_mode: shortest`
- `path_mode: chemical_shortest`
- `excluded_species_formulas`
- `excluded_species_labels`
- `hub_penalty_strength`
- `max_paths`

The `chemical_shortest` mode is especially useful when highly connected public intermediates such as `H` would otherwise dominate shortest-path searches.

## 11. Typical Use Cases

### 11.1 Metal surface catalysis

In standard slab-based catalysis, the most important pieces are often:

- host definition
- geometric-site annotation
- normal species/reaction analysis

### 11.2 Zeolite acid-center systems

In zeolite and acid-center systems, the most important layer is often `active_site`, because the key scientific questions usually concern:

- BAS/LAS identity
- state changes of acid centers
- coupling between molecular reactions and site-state evolution

### 11.3 Local metal centers in hosts

Systems containing local centers such as `Ga(OH)2` are particularly suitable for the active-site framework.

## 12. Practical Usage Pattern

A practical workflow is often:

1. Run the main analyzer and generate the full bipartite graph.
2. Inspect `frame_species.csv`, transform events, and active-site outputs.
3. Use `focus` to obtain a readable backbone graph.
4. Use `context` for local inspection around a known node.
5. Use `story` for source-to-target route extraction.

## 13. Current Limitations

Important current limits include:

1. Site-state labels are still rule-based, not a full chemical naming engine.
2. `active_site` definitions still require explicit chemical judgment from the user.
3. `story` now supports `shortest` and `chemical_shortest`, but it still does not provide more advanced route-family enumeration or richer path-scoring schemes.
4. Graph postprocessing is intentionally centered on the bipartite graph.
5. Some rules are still specialized to particular chemical families such as BAS and metal-centered LAS cases.

## 14. Development Suggestions

Promising future directions include:

1. More chemically expressive site-state labels.
2. Stronger route-ranking strategies for `story`.
3. Higher-level automated summaries and reports.
4. Better semi-automatic active-site definition tools.
5. Minimal smoke tests to protect future refactors.

## 15. Summary

`rxn_analyzer` has evolved into a structured catalytic-trajectory analysis framework rather than a simple CSV-export script. Its main strengths are:

- converting trajectories into reaction events
- supporting both normal species analysis and site-participation analysis
- distinguishing static geometric sites from dynamic active sites
- exporting layered outputs suitable for scripting, manual reading, and graph visualization
- providing task-oriented graph postprocessing for backbone extraction, local context extraction, and route-story extraction
