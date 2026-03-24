# rxn_analyzer 用户手册

## 1. 文档目的

本文档用于系统说明 `rxn_analyzer` 程序的设计目标、总体架构、核心对象、主要算法思路、使用方法、输出结果解释、典型适用场景与后处理分析方式。文档面向以下几类读者：

- 希望直接运行程序并分析轨迹结果的研究人员
- 需要理解程序设计思路与边界条件的开发者
- 需要评估程序是否适用于某类催化体系的课题负责人

`rxn_analyzer` 的核心定位不是一个通用分子动力学后处理工具箱，而是一个围绕“反应网络提取”和“位点参与分析”设计的专门分析程序。它尤其适合如下任务：

- 从轨迹中识别成键与断键事件
- 对非宿主物种进行逐帧识别与跨帧映射
- 构建反应事件表与反应网络图
- 在静态几何位点和动态活性位点两个层面分析反应发生位置
- 研究宿主骨架上的活性中心是否直接参与反应

在当前版本中，程序已经支持将分析对象分为三个概念层次：

- `host`：宿主骨架或基底背景
- `geometric_site`：静态几何位点
- `active_site`：动态活性中心

这三个层次共同构成了程序的概念框架。

---

## 2. 程序设计目标与总体思路

### 2.1 为什么需要这个程序

对于反应轨迹，尤其是催化反应轨迹，最常见的困难不是“有没有原子坐标”，而是“如何把连续轨迹变成可解释的化学事件”。一个轨迹文件本身只告诉我们每一帧有哪些原子、坐标是多少，但并不会自动回答以下问题：

- 哪些键形成了，哪些键断裂了
- 某个中间体是在什么时候出现的
- 一个分子是气相态、吸附态还是位点本体的一部分
- 某个反应发生在什么位点上
- 位点是否只是提供环境，还是本身也发生了化学变化

`rxn_analyzer` 的设计目标，就是把“轨迹坐标数据”逐步转化为“可用于化学讨论的事件和网络”。

### 2.2 程序的核心思想

程序的核心思想可以概括为三步：

1. 从几何轨迹恢复化学连接关系  
   程序通过阈值、滞后策略和持久化判定，识别每一帧中的共价边、吸附边与宿主边，并记录真正可信的成键/断键事件。

2. 从连接关系识别化学对象  
   程序把非 `host` 原子组成的连通分量识别为 `species`，并为每个组分赋予标签，如分子式、SMILES、吸附签名、几何位点签名等。

3. 从相邻帧变化提取反应语义  
   程序将相邻帧中的物种变化、键变化和位点状态变化综合起来，提取普通反应事件、活性位点事件、联合反应事件以及图结构化结果。

### 2.3 为什么要区分 geometric_site 和 active_site

这是程序设计里非常关键的一点。

在很多旧式分析框架里，“位点”只有一个概念，既可以指金属表面的 top、bridge、fcc，也可以指分子筛中的酸中心。但在化学上，这两种“位点”其实并不属于同一类对象：

- `geometric_site` 更像“空间位置标签”
- `active_site` 更像“化学活性中心对象”

因此本程序明确把两者分开：

- `geometric_site` 主要回答：这个分子挂在哪个几何位置上
- `active_site` 主要回答：这个活性中心当前处于什么状态，它是否参与了反应

这种区分让程序同时适用于：

- 金属表面催化体系
- 分子筛 / 酸中心体系
- 含有局部金属活性中心的宿主-客体体系

---

## 3. 核心概念与对象模型

### 3.1 host

`host` 表示宿主骨架、基底背景、固定框架或大尺度支撑环境。  
在许多体系中，`host` 对应：

- 金属 slab
- 分子筛骨架
- 固定框架氧化物
- 大尺寸宿主结构

程序中 `host` 的作用主要有三个：

1. 定义哪些原子不属于普通自由物种
2. 定义吸附边和宿主边的类别
3. 为 `geometric_site` 与 `active_site` 提供锚点

当前 `host_definition` 支持：

- 按元素定义
- 按索引定义
- 索引反选
- 区间简写
- `0-based / 1-based` 索引方式

这使得程序既适合简单 slab，也适合复杂、局部定义的宿主体系。

### 3.2 species

`species` 是程序中对非 `host` 连通分量的抽象。  
它代表的是：

- 气相分子
- 吸附中间体
- 产物
- 反应过程中出现的任何非宿主连通组分

每个 `species` 都会在单帧中被赋予标签，标签来源包括：

- 分子式
- SMILES 或 WL 哈希
- 吸附签名
- 几何位点签名

程序通过跨帧映射把当前帧与上一帧的 `species` 对齐，从而判断是否发生了反应。

### 3.3 geometric_site

`geometric_site` 是静态定义的位点模板，通常由用户通过一组 `host` 原子来定义。  
它适合：

- top
- bridge
- hollow
- 四配位位点
- 某些固定孔口位

`geometric_site` 的本质是：  
给 `species` 提供一个“这个物种当前挂在哪个位置”的附加标签。

它本身默认不演化，也不承担“位点本体状态变化”的分析任务。

### 3.4 active_site

`active_site` 是程序中最有特色的一层。  
它不是一个纯几何位置，而是一个可演化的化学中心。

每个 `active_site` 由两层构成：

- `core_members`
  由 `host` 原子组成，用来锚定“这个位点是谁”
- `state_members`
  描述位点当前状态的可变成员，包括：
  - `intrinsic`
  - `incorporated`
  - `associated`

这三种成员的含义是：

- `intrinsic`
  本来就属于位点本体的部分
- `incorporated`
  原本可能是外来片段，但现在已经并入位点状态
- `associated`
  与位点强关联，但尚未并入位点本体

这种设计使程序能表达：

- 位点只是吸附了某个分子
- 位点吸收了一个新的局部片段
- 位点本身发生了化学状态变化

---

## 4. 程序总体架构

### 4.1 架构总览

程序目前可以概括为五层核心运行模块：

1. `config_loader`
2. `analyzer`
3. `species`
4. `active_site`
5. `graph`

外加若干辅助输出与后处理模块。

### 4.2 config_loader

`config_loader` 负责读取 YAML 配置、校验字段、构建运行时对象。  
当前版本只接受新命名配置：

- `host_definition`
- `geometric_site`
- `active_site`

旧命名如 `slab_definition`、`site`、`reactive_site` 已不再兼容。

### 4.3 analyzer

`analyzer` 是总调度层。  
它自身尽量不承载复杂业务细节，而是负责协调以下流程：

1. 当前帧进入分析器
2. `EdgePipeline` 识别当前连接关系与键事件
3. `SpeciesPipeline` 生成当前 `SpeciesFrameSnapshot`
4. `SpeciesRuntime` 进行跨帧映射与普通反应事件发射
5. `ActiveSiteRuntime` 识别位点状态并发射位点事件
6. `OutputWriter` / `ActiveSiteOutputWriter` 写出最终输出

### 4.4 species 子模块

`species/` 已被拆成独立子模块，主要包括：

- `chem.py`
  单帧物种识别相关的基础化学工具
- `model.py`
  `SpeciesFrameSnapshot`
- `pipeline.py`
  单帧物种识别
- `mapping.py`
  跨帧 component 对齐
- `emitter.py`
  普通 transform 事件发射
- `runtime.py`
  species 主运行时串联

这条线主要处理“普通非宿主反应”。

### 4.5 active_site 子模块

`active_site/` 是另一条完整的独立分析线，主要包括：

- `model.py`
  活性位点核心数据模型
- `types.py`
  运行时中间对象
- `pipeline.py`
  单帧位点状态识别
- `tracker.py`
  跨帧位点事件跟踪
- `coupling.py`
  位点状态变化与普通反应的耦合
- `runtime.py`
  active-site 主运行时
- `output.py`
  active-site 输出与 site-aware graph
- `rules.py`
  位点状态标签规则

这条线主要处理“位点自身状态演化”和“位点与反应的关系”。

### 4.6 graph 子模块

`graph/` 负责反应网络构建与 GraphML 属性管理。  
目前核心包括：

- `ids.py`
  节点与反应 ID 分配
- `attrs.py`
  节点属性与标签字段解析
- `network.py`
  主二部反应图构建
- `postprocess/`
  图后处理模式

当前版本程序只输出主二部图，不再输出旧式普通图。

---

## 5. 逐帧分析流程

### 5.1 边与键事件识别

每一帧首先进入 `EdgePipeline`。  
这一层会识别：

- 共价边
- 吸附边
- 宿主-宿主边
- 成键与断键事件

程序使用阈值、滞后和持久化策略，尽量避免把热波动误当作真正化学反应。

### 5.2 species 快照构建

接着 `SpeciesPipeline` 会把所有非 `host` 原子构成的连通分量识别出来，并为其赋予标签。  
结果会整理成一个 `SpeciesFrameSnapshot`，其中包含：

- components
- labels
- multiset
- component labels
- 吸附配对信息
- 几何位点赋值

### 5.3 普通反应事件提取

`SpeciesRuntime` 会将当前快照与上一帧比较，结合 component mapping，提取：

- 物种消失
- 物种出现
- 分裂/合并
- 普通 transform events

这些事件最终进入：

- `events_transform.csv`
- 主二部图

### 5.4 active site 状态构建

与 species 并行，`ActiveSiteRuntime` 会分析当前帧的每个 `active_site`。  
它会完成：

- 核心成员校验
- intrinsic/incorporated/associated 分类
- 状态组分构造
- `state_formula`
- `state_topology`
- `state_label`

当前程序已经支持一定程度的化学语义标签，例如：

- `ga_h`
- `ga_h2`
- `ga_oh`
- `ga_oh2`
- `ga_h_oh`
- `protonated_bas`
- `deprotonated_bas`
- `alkoxy_bas`

### 5.5 active site 事件与耦合

在跨帧层面，程序会进一步识别：

- `site_state_change`
- `joint_site_reactions`
- `site_reaction_coupling`

由此可以区分三类重要情况：

- 反应发生，但位点状态不变
- 反应发生，同时位点状态变化
- 没有明显普通反应，但位点自身状态发生变化

这正是活性中心分析的核心。

---

## 6. 配置文件说明

### 6.1 主分析配置

最常用的主配置项包括：

- `run`
  输入轨迹、输出前缀、步长、最大帧数、进度条设置
- `criteria`
  成键/断键阈值与持久化参数
- `host_definition`
  宿主定义
- `geometric_site`
  静态位点配置
- `active_site`
  活性位点配置
- `analyzer`
  分析器行为参数

### 6.2 host_definition

可以按元素：

```yaml
host_definition:
  mode: "elements"
  elements: ["Pt"]
```

也可以按索引：

```yaml
host_definition:
  mode: "indices"
  indices: ["1-20", 22, "26-100"]
  index_base: 1
  invert: false
```

如果希望“没有 host”，也可以：

```yaml
host_definition:
  mode: "indices"
  indices: []
  index_base: 0
  invert: false
```

### 6.3 geometric_site

几何位点通过独立 YAML 文件定义。  
用户可在 `analyzer.geometric_site_signature_mode` 中决定是否把位点签名写进 species 标签。

常见模式：

- `none`
- `type`
- `id`
- `type+id`

### 6.4 active_site

活性位点同样通过独立 YAML 文件定义。  
典型示例：

```yaml
- id: "B_Ga_01"
  family: "BAS_Ga"
  core_members: [616, 31, 35]
  initial_state_members: [608]
  allowed_state_elements: ["H", "O", "C"]
  max_state_members: 20
  rule_profile: "bas_default"
```

对于分子筛酸中心体系，这种定义方式非常灵活，因为它允许：

- 用 `core_members` 锚定位点身份
- 用 `initial_state_members` 定义初始本体状态
- 再由程序动态识别后续状态变化

---

## 7. 主要输出文件

### 7.1 普通事件与主图

常见主输出包括：

- `<out_prefix>_events_transform.csv`
- `<out_prefix>_events_bonds.csv`
- `<out_prefix>_network.graphml`
- `<out_prefix>_reactions_summary.tsv`
- `<out_prefix>_frame_species.csv`

其中：

- `events_transform.csv`
  是普通反应事件表
- `events_bonds.csv`
  是成键/断键证据表
- `network.graphml`
  是主二部反应图
- `reactions_summary.tsv`
  是正逆反应汇总表
- `frame_species.csv`
  是逐帧物种快照

### 7.2 active_site 输出

启用 `active_site` 后，还会增加：

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

这些输出分别承担：

- 状态轨迹
- 状态变化事件
- 联合反应语义
- 普通反应与位点状态耦合
- 状态统计 summary
- family-level 聚合 summary
- site-aware graph 可视化

---

## 8. 图结构与后处理

### 8.1 主图结构

当前主图采用二部图：

- `species` 节点
- `reaction` 节点

边表示：

- 物种作为反应物参与某个反应
- 某个反应产生某个物种

这种结构比旧式普通图更适合表达化学语义，也更适合后续做正逆合并、反应筛选和路径提取。

### 8.2 site-aware graph

除了主图外，程序还可以输出 `site-aware graph`，其节点类型包括：

- `species`
- `reaction`
- `site_state`

它用于表达：

- 某个反应是否发生在某个活性位点上
- 该反应是否改变了位点状态

### 8.3 graph postprocess

图后处理现在已经重构为三种任务模式：

- `focus`
- `context`
- `story`

所有模式都使用 YAML 配置。

#### focus

`focus` 用于提取主干图。  
适合回答：

- 网络中最重要的节点和反应是什么

当前支持：

- `score_mode`
  - `weighted_degree`
  - `reaction_weight`
  - `hybrid`
- `top_species`
- `top_reactions`
- `top_species_per_family`
- `species_family_mode`
- `min_edge_weight`
- `collapse_reversible`
- `include_neighbors`
- `keep_nodes`
- `protect_seed_neighbors`

其中：

- `top_species_per_family + species_family_mode=ads_state`
  可以在 `ads` 和 `non_ads` 两类之间保留更多平衡
- `protect_seed_neighbors`
  用于围绕手动指定的种子节点保留局部上下文

#### context

`context` 用于围绕一个或多个 seed 提取局部子图。  
适合回答：

- 围绕某个中间体，附近发生了什么

支持：

- `seeds`
- `depth`
- `direction`
- `min_edge_weight`
- `collapse_reversible`
- `prune_isolates`

#### story

`story` 用于提取 source 到 target 的路径子图。  
适合回答：

- 从反应物 A 到产物 B 的主要路径是什么

当前第一版支持：

- `sources`
- `targets`
- `path_mode: shortest`
- `direction`
- `max_paths`
- `min_edge_weight`
- `collapse_reversible`
- `prune_isolates`

---

## 9. 适用场景

### 9.1 金属表面催化

在金属表面体系中，最重要的通常是：

- 哪些中间体出现
- 哪些反应最常见
- 中间体吸附在哪类几何位点

这时建议主要启用：

- `host`
- `geometric_site`

而 `active_site` 通常不是重点。

### 9.2 分子筛酸中心体系

在分子筛或局部酸中心体系中，更关心的是：

- 某个 BAS / LAS 当前处于什么状态
- 位点是否参与反应
- 位点状态是否变化

这时建议主要启用：

- `host`
- `active_site`

`geometric_site` 往往不是主工具。

### 9.3 含局部金属中心的宿主体系

例如：

- `Ga(OH)2`
- 单金属位
- 框架外金属酸中心

这类体系非常适合 `active_site` 分析，因为程序已经支持：

- `ga_h`
- `ga_h2`
- `ga_oh`
- `ga_oh2`
- `ga_h_oh`

等更化学化的状态标签。

---

## 10. 典型使用流程

### 10.1 仅做普通反应分析

关闭：

- `geometric_site`
- `active_site`

得到：

- 普通物种事件
- 主反应图

### 10.2 普通反应 + 几何位点标注

开启：

- `geometric_site`

关闭：

- `active_site`

得到：

- 物种标签中包含几何位点信息

### 10.3 活性位点状态分析

关闭：

- `geometric_site`

开启：

- `active_site`

得到：

- 位点状态表
- 位点事件表
- site-aware graph

### 10.4 图后处理

对生成出的 `network.graphml` 再使用：

- `focus`
  提主干
- `context`
  看局部
- `story`
  看路径

这使程序从“生成原始网络”进一步扩展到“生成可读、可汇报的简化网络”。

---

## 11. 局限性与当前边界

虽然程序已经具备较完整的分析链，但当前仍有一些边界需要注意：

1. 位点状态标签仍是规则驱动  
   程序可以给出较好的工程化标签，但并不是自动量子化学命名器。

2. `active_site` 需要用户显式定义  
   虽然已有辅助生成脚本，但位点定义本身仍需要化学判断。

3. 后处理 `story` 目前只支持最短路径  
   还没有加入更复杂的路径评分模式。

4. 图后处理仍主要围绕二部图  
   当前设计是明确偏向二部图的，这和程序主输出是一致的。

5. 某些规则高度依赖体系类型  
   例如 `BAS`、`LAS_Ga` 的状态规则是为特定化学场景优化的，后续仍可能需要继续扩展。

---

## 12. 后续发展建议

对于后续版本，我建议优先考虑以下方向：

1. 增强位点状态标签的化学语义
2. 为 `story` 模式加入高权重路径或 top-k 路径策略
3. 增加更高层的自动 summary / report 输出
4. 继续补充 `active_site` 自动定义工具
5. 增加最小 smoke tests 以保护后续重构

---

## 13. 总结

`rxn_analyzer` 当前已经不只是一个“从轨迹中生成几个 CSV”的脚本，而是一个较完整的催化反应轨迹分析框架。它的核心优势在于：

- 将轨迹几何信息转化为反应事件
- 同时支持普通物种分析与位点参与分析
- 区分静态几何位点与动态活性位点
- 输出适合脚本处理、人工阅读和图可视化的多层结果
- 提供了后处理机制，将复杂网络进一步收缩成主干、局部上下文和路径故事

对于金属表面催化、分子筛酸中心体系以及带有局部活性中心的宿主体系，它都已经具备较好的实用性。随着规则和 summary 的继续增强，这个程序可以进一步发展为一个更加稳定、更加化学友好的反应网络分析平台。
