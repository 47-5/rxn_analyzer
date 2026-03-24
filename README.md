# rxn_analyzer

`rxn_analyzer` 用于从轨迹中识别成键/断键、构建物种标签、提取反应事件，并进一步分析宿主骨架、几何位点和活性位点。

当前推荐使用的新命名是：

- `host`
  表示宿主骨架或基底背景。历史名称是 `slab`。
- `geometric_site`
  表示静态几何位点，例如 top / bridge / fcc。历史名称是 `site`。
- `active_site`
  表示会参与反应、且状态可能随时间变化的活性中心。历史名称是 `reactive_site`。

旧命名不再兼容。配置、CLI 和文档现在统一只使用新名字。

## 示例文件

推荐先看这几个文件：

- [example_analyzer.yaml](/D:/code/rxn_analyzer/example_config/example_analyzer.yaml)
- [example_analyzer_clean.yaml](/D:/code/rxn_analyzer/example_config/example_analyzer_clean.yaml)
- [example_site.yaml](/D:/code/rxn_analyzer/example_config/example_site.yaml)
- [example_reactive_site.yaml](/D:/code/rxn_analyzer/example_config/example_reactive_site.yaml)

其中：

- `example_analyzer.yaml`
  完整示例配置，注释最详细。
- `example_analyzer_clean.yaml`
  更适合直接复制修改的精简版，已切到新命名。
- `example_site.yaml`
  静态几何位点定义示例，对应 `geometric_site`。
- `example_reactive_site.yaml`
  动态活性位点定义示例，对应 `active_site`。

## 核心对象

### 1. host

`host` 是宿主骨架定义，对应配置段：

```yaml
host_definition:
  mode: "elements"
  elements: ["Pt"]
```

它定义哪些原子属于宿主背景。

如果你想按索引定义 `host`，现在还支持：

- 反选模式 `invert: true`
- 区间简写，例如 `["1-20", 22, "26-100"]`
- 选择 `index_base: 0` 或 `index_base: 1`

例如：

```yaml
host_definition:
  mode: "indices"
  indices: ["1-20", 22, "26-100"]
  index_base: 1
  invert: false
```

### 2. geometric_site

`geometric_site` 适合静态几何位点：

- top
- bridge
- fcc
- hcp
- fourfold

它的核心作用是给普通物种附加“吸附在哪个几何位点上”的标签。

对应配置段：

```yaml
geometric_site:
  enabled: true
  file: "example_site.yaml"

analyzer:
  geometric_site_signature_mode: "type"
```

### 3. active_site

`active_site` 适合会参与反应、并且状态会变化的活性中心：

- 分子筛 BAS / LAS
- 单原子位点
- 团簇位点
- 缺陷位点
- 像 `Ga(OH)2` 这类局部化学中心

对应配置段：

```yaml
active_site:
  enabled: true
  file: "example_reactive_site.yaml"

analyzer:
  active_site_record_states: true
  active_site_record_events: true
  active_site_record_joint_reactions: true
```

## 什么时候用 geometric_site，什么时候用 active_site

### 金属表面 / 标准 slab 体系

通常更关心：

- 分子挂在哪个 top / bridge / hollow 上
- 不太关心位点自身化学身份变化

这时优先用：

- `host`
- `geometric_site`

### 分子筛 / 酸中心体系

通常更关心：

- 某个 BAS / LAS 当前是什么状态
- 位点是否参与反应
- 反应是否改变了位点状态

这时优先用：

- `host`
- `active_site`

`geometric_site` 往往不是主工具。

## 推荐使用路径

### 1. 只做普通反应网络分析

```yaml
geometric_site:
  enabled: false

active_site:
  enabled: false

analyzer:
  geometric_site_signature_mode: "none"
```

### 2. 做普通反应分析 + 几何位点标注

```yaml
geometric_site:
  enabled: true
  file: "example_site.yaml"

active_site:
  enabled: false

analyzer:
  geometric_site_signature_mode: "type+id"
```

### 3. 做活性位点状态分析

```yaml
geometric_site:
  enabled: false

active_site:
  enabled: true
  file: "example_reactive_site.yaml"
```

### 4. 同时保留几何位点和活性位点

```yaml
geometric_site:
  enabled: true
  file: "example_site.yaml"

active_site:
  enabled: true
  file: "example_reactive_site.yaml"
```

## 运行方式

推荐直接使用模块入口：

```bash
python -m rxn_analyzer.runner --config example_config/example_analyzer_clean.yaml
```

Windows 上也可以：

```bash
py -3 -m rxn_analyzer.runner --config example_config/example_analyzer_clean.yaml
```

## 代码结构

当前代码里最值得关注的两个子模块是：

- [`rxn_analyzer/species/`](/D:/code/rxn_analyzer/rxn_analyzer/species)
  负责普通非 host 物种识别与反应解释
- [`rxn_analyzer/active_site/`](/D:/code/rxn_analyzer/rxn_analyzer/active_site)
  负责动态活性位点状态、事件和耦合分析

其中：

- [`rxn_analyzer/species/chem.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/chem.py)
  单帧 chemistry 工具，例如连通分量、formula、SMILES、ads 签名
- [`rxn_analyzer/species/model.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/model.py)
  `SpeciesFrameSnapshot`
- [`rxn_analyzer/species/pipeline.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/pipeline.py)
  单帧物种识别
- [`rxn_analyzer/species/mapping.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/mapping.py)
  跨帧 component 对齐
- [`rxn_analyzer/species/emitter.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/emitter.py)
  transform 事件发射
- [`rxn_analyzer/species/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/runtime.py)
  species 主运行时串联

- [`rxn_analyzer/active_site/model.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/model.py)
  `ActiveSite` 相关核心模型
- [`rxn_analyzer/active_site/types.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/types.py)
  active-site 运行时中间类型
- [`rxn_analyzer/active_site/pipeline.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/pipeline.py)
  单帧 active-site 状态识别
- [`rxn_analyzer/active_site/tracker.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/tracker.py)
  跨帧 active-site 事件生成
- [`rxn_analyzer/active_site/coupling.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/coupling.py)
  active-site 与普通反应的耦合分析
- [`rxn_analyzer/active_site/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/runtime.py)
  active-site 主运行时串联

主调度仍在：

- [`rxn_analyzer/analyzer.py`](/D:/code/rxn_analyzer/rxn_analyzer/analyzer.py)
- [`rxn_analyzer/runner.py`](/D:/code/rxn_analyzer/rxn_analyzer/runner.py)

如果你要继续开发，推荐的阅读顺序是：

1. [`rxn_analyzer/analyzer.py`](/D:/code/rxn_analyzer/rxn_analyzer/analyzer.py)
2. [`rxn_analyzer/species/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/species/runtime.py)
3. [`rxn_analyzer/active_site/runtime.py`](/D:/code/rxn_analyzer/rxn_analyzer/active_site/runtime.py)
4. 再按需要深入 `pipeline / mapping / emitter / tracker / coupling`

更完整的开发向说明见：

- [ARCHITECTURE.md](/D:/code/rxn_analyzer/ARCHITECTURE.md)

## 自动生成 active_site 定义

如果你在分子筛体系里有很多类似的 BAS/LAS，不想手动一个个写 `active_site`，现在可以先用自动生成工具从首帧结构里提取。

当前第一版优先支持你这种 `center - O(H) - framework` 的 BAS 模式，例如：

- `Ga - O(H) - Si`

辅助脚本：

```bash
py -3 generate_active_sites.py
```

这个脚本不是包内正式 API，而是一个独立的启发式辅助工具。推荐直接在 IDE 里打开 [`generate_active_sites.py`](/D:/code/rxn_analyzer/generate_active_sites.py)，修改 `if __name__ == "__main__":` 下面的参数块再运行。

现在推荐的用法是：

- `CONFIG_PATH`
  只用来读取 `host_definition` 和 `criteria`
- `STRUCTURE_PATH`
  直接指向首帧结构文件，优先使用这个文件，而不是整条轨迹
- `TRAJ_OVERRIDE`
  只有在 `STRUCTURE_PATH = None` 时才作为后备方案使用

也就是说，如果你只是想快速生成 `active_site` 种子，最省时的做法是保留原分析配置，再单独提供一个首帧结构文件给 `STRUCTURE_PATH`。

例如你可以这样写：

```python
CONFIG_PATH = "GaH2_gfn_xtb_analyzer.yaml"
STRUCTURE_PATH = "GaH2_frame0.xyz"
CENTER_ELEMENTS = ["Ga", "Al"]
FRAMEWORK_ELEMENTS = ["Si"]
ID_PREFIX = "B"
FAMILY = "BAS"
PER_CENTER_NAMING = True
```

这样会自动生成类似：

- `B_Ga_01`, `B_Ga_02`, ...
- `B_Al_01`, `B_Al_02`, ...

对应的 `family` 也会自动变成：

- `BAS_Ga`
- `BAS_Al`

这会自动寻找类似下面这种定义：

```yaml
- id: "B_Ga_01"
  family: "BAS_Ga"
  core_members: [Ga, O, Si]
  initial_state_members: [H]
```

目前这个工具仍然是“半自动”的：

- 它非常适合先帮你批量起草 active site
- 生成后仍建议你快速人工抽查几条
- 后面如果你还需要别的 motif，可以继续往这个工具里加

## 常见输出

主输出通常包括：

- `<out_prefix>_events_transform.csv`
- `<out_prefix>_events_bonds.csv`
- `<out_prefix>_network.graphml`
- `<out_prefix>_reactions_summary.tsv`
- `<out_prefix>_frame_species.csv`

启用 `active_site` 后，还会增加：

- `<out_prefix>_site_states.csv`
- `<out_prefix>_site_events.csv`
- `<out_prefix>_joint_site_reactions.csv`
- `<out_prefix>_site_reaction_coupling.csv`
- `<out_prefix>_network_site_aware.graphml`
- `<out_prefix>_site_outputs_README.txt`

注意：

- 这些输出文件名目前仍保留历史上的 `site` / `reactive_site` 风格命名。
- 配置和概念层面推荐使用 `host / geometric_site / active_site`。

## `frame_species.csv` 里和 active_site 相关的列

如果开启了 `record_frame_species`，你会看到这些和活性位点有关的列：

- `active_site_roles`
  每个连通分量在活性位点视角下的主角色，可能是 `intrinsic`、`incorporated`、`associated`、`mixed` 或 `none`。
- `active_site_ids`
  该分量关联到哪些活性位点。
- `active_site_summaries`
  适合直接阅读的简短摘要，例如 `LAS_Ga_01:intrinsic`。
- `is_active_site_owned`
  是否属于某个活性位点本体，也就是 `intrinsic/incorporated`。
- `is_active_site_associated`
  是否只是与活性位点相关联的外来分量。
- `active_site_memberships`
  详细 JSON 明细，用于精确追踪原子级归属。

## 关于旧命名

以下旧名字已经废弃，不能再用于新版本配置：

- `slab_definition`
- `site`
- `reactive_site`
- `analyzer.site_signature_mode`
- `analyzer.reactive_site_record_*`

请统一改用：

- `host_definition`
- `geometric_site`
- `active_site`
- `analyzer.geometric_site_signature_mode`
- `analyzer.active_site_record_*`
- `analyzer.active_site_streaming`

## Graph Postprocess

Graph postprocessing now lives under:

- `rxn_analyzer.graph.postprocess`

The new interface is YAML-driven and task-oriented. The first supported mode is:

- `focus`
  Keep the most important `species` and `reaction` nodes from an existing GraphML network.
- `context`
  Extract a local neighborhood around one or more seed nodes.
- `story`
  Extract a source-to-target path subgraph.

Example config:

- [example_graph_postprocess_focus.yaml](/D:/code/rxn_analyzer/example_config/example_graph_postprocess_focus.yaml)
- [example_graph_postprocess_context.yaml](/D:/code/rxn_analyzer/example_config/example_graph_postprocess_context.yaml)
- [example_graph_postprocess_story.yaml](/D:/code/rxn_analyzer/example_config/example_graph_postprocess_story.yaml)

Key `focus` options:

- `score_mode`
  Controls how node importance is scored before keeping the top nodes.
  Supported values:
  - `weighted_degree`
    Rank nodes by the sum of incident edge weights.
  - `reaction_weight`
    Rank reaction nodes by their own reaction weight, and rank species by the total weight of neighboring reaction nodes.
    This is usually the most chemical-intuitive option for the current bipartite graph.
  - `hybrid`
    Sum `weighted_degree` and `reaction_weight`.
- `top_species`
  Number of species nodes kept before optional neighbor expansion.
- `top_reactions`
  Number of reaction nodes kept before optional neighbor expansion.
- `top_species_per_family`
  Optional quota-style protection against one species family taking over all retained species slots.
  First keep up to this many top-ranked species from each family, then fill the remaining `top_species` slots globally.
  Set to `0` to disable.
- `species_family_mode`
  Current supported value:
  - `ads_state`
    Split species into `ads` and `non_ads`.
- `min_edge_weight`
  Remove bipartite edges whose `weight` is below this threshold before scoring.
- `collapse_reversible`
  For the bipartite graph, merge true forward/reverse reaction-node pairs into one reversible reaction node.
  Single-direction repeated reactions stay directional.
- `include_neighbors`
  After selecting the top-ranked nodes, also keep their first neighbors to preserve local context.
- `keep_nodes`
  Extra node labels / `orig_id` values that must be kept even if their score is low.
- `protect_seed_neighbors`
  Protect the nodes listed in `keep_nodes` together with their N-hop neighborhood, even if their score is low.

Key `context` options:

- `seeds`
  Seed nodes to expand from. Each seed can be matched by node id, `orig_id`, `label`, or `display_label`.
- `depth`
  Number of graph hops to expand from the seeds.
- `direction`
  One of:
  - `both`
  - `out`
  - `in`
- `min_edge_weight`
  Remove weak edges before context expansion.
- `collapse_reversible`
  Optional preprocessing step before extracting the local neighborhood.
- `prune_isolates`
  Remove isolated nodes after extracting the context subgraph.

Key `story` options:

- `sources`
  Source nodes for path extraction.
- `targets`
  Target nodes for path extraction.
- `path_mode`
  Current supported values:
  - `shortest`
  - `chemical_shortest`
    Use weighted shortest paths that penalize hub-like intermediate species.
- `direction`
  One of:
  - `directed`
  - `undirected`
- `max_paths`
  Maximum number of shortest paths retained in the story subgraph.
- `hub_penalty_strength`
  Penalty strength used by `chemical_shortest` to downweight highly connected intermediate species.
- `excluded_species_formulas`
  Explicitly remove matching species formulas from the path search, except when they are sources or targets.
- `excluded_species_labels`
  Explicitly remove matching species labels from the path search, except when they are sources or targets.
- `min_edge_weight`
  Remove weak edges before path extraction.
- `collapse_reversible`
  Optional preprocessing step before path extraction.
- `prune_isolates`
  Remove isolated nodes after extracting the story subgraph.

Recommended usage:

```bash
python -m rxn_analyzer.graph.postprocess.runner --config example_config/example_graph_postprocess_focus.yaml
```

If you prefer the same style as `run_analyze.py`, use the root script:

```bash
python run_postprocess.py
```
