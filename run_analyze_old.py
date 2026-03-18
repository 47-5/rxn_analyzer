from tqdm import tqdm
from rxn_analyzer import (
    Criteria, DistanceHysteresisParams,
    SlabDefinition, ReactionAnalyzer, AnalyzerConfig
)
from rxn_analyzer.io import frames
from rxn_analyzer import network as net

net.reset_node_id_mapping()


def count_frames(path: str, stride: int = 1) -> int:
    # 用你现成的 frames()，计数逻辑与实际处理完全一致
    return sum(1 for _k, _atoms in frames(path, stride=stride))


if __name__ == "__main__":
    traj = "C12H26.extxyz"
    stride = 1

    criteria = Criteria(
        persist=5,
        cooldown=3,
        cov=DistanceHysteresisParams(1.15, 1.25, None),
        ads=DistanceHysteresisParams(1.25, 1.40, None),
        slab=DistanceHysteresisParams(1.35, 1.55, 12),
    )

    slab = SlabDefinition(elements=frozenset({"Pt"}))

    config = AnalyzerConfig(
        wl_iters=3,
        baseline_persist=5,
        record_init_events=False,
        drop_init_events=True,
        ads_signature_mode="coarse",

        smiles_recompute_mode="on_change",

        smiles_mode="rdkit_topology",
        smiles_combo_treat_suspicious_as_failure=True,

        rdkit_sanitize=True,
        rdkit_hide_hs=False,
        smiles_fallback_to_formula_if_suspicious=False,

        # --- component overlap split ---
        split_by_component=True,
        component_overlap_min=0.0,
        component_overlap_mode="symmetric",

        record_frame_species=True,
        frame_species_streaming=True,

        # --- 控制 *_reactions_summary.tsv 是否包含正/逆反应帧列表 ---
        reaction_summary_include_frames=True,
    )

    an = ReactionAnalyzer(criteria=criteria, slab_def=slab, config=config)

    total = count_frames(traj, stride=stride)

    for frame, (_k, atoms) in enumerate(
        tqdm(frames(traj, stride=stride), total=total, desc="Analyzing", unit="frame")
    ):
        an.process_frame(frame, atoms)

    an.write_outputs("C12H26")