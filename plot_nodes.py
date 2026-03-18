#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import xml.etree.ElementTree as ET
from io import BytesIO

from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


SMILES_RE = re.compile(r"smiles=([^|]+)")


def load_graphml_nodes(path: Path):
    """
    Return:
      nodes: list of (node_id, data_map[key_id] = text)
      key_id_to_name, key_name_to_id
    """
    tree = ET.parse(path)
    root = tree.getroot()
    ns = {"g": "http://graphml.graphdrawing.org/xmlns"}

    key_id_to_name: Dict[str, str] = {}
    key_name_to_id: Dict[str, str] = {}

    for key in root.findall("g:key", ns):
        kid = key.get("id")
        name = key.get("attr.name")
        if kid and name:
            key_id_to_name[kid] = name
            key_name_to_id[name] = kid

    graph = root.find("g:graph", ns)
    if graph is None:
        raise ValueError("GraphML 中未找到 <graph> 元素")

    nodes = []
    for node in graph.findall("g:node", ns):
        nid = node.get("id", "")
        data_map: Dict[str, str] = {}
        for data in node.findall("g:data", ns):
            k = data.get("key")
            if not k:
                continue
            data_map[k] = (data.text or "")
        nodes.append((nid, data_map))

    return nodes, key_id_to_name, key_name_to_id


def resolve_key(user_key: str, key_id_to_name: Dict[str, str], key_name_to_id: Dict[str, str]) -> str:
    """
    user_key can be key id like 'd13' or attr.name like 'label_smiles'
    returns key id
    """
    if user_key in key_id_to_name:
        return user_key
    if user_key in key_name_to_id:
        return key_name_to_id[user_key]
    raise KeyError(f"未找到 key: {user_key}（既不是 key id 也不是 attr.name）")


def extract_smiles_candidates(text: str, delimiter: str = " + ") -> List[str]:
    """
    从 data 文本中提取候选 SMILES。
    - 若出现 smiles=...，只使用 smiles= 的值（更可靠）
    - 否则回退到普通片段
    - 若包含 | ，只保留 | 左侧
    """
    if not text:
        return []

    t = " ".join(text.split())

    # 优先用 smiles=
    if "smiles=" in t:
        out = []
        for m in SMILES_RE.findall(t):
            m = m.strip()
            if "|" in m:
                m = m.split("|", 1)[0].strip()
            if m:
                out.append(m)
        # 去重保序
        seen = set()
        uniq = []
        for s in out:
            if s not in seen:
                seen.add(s)
                uniq.append(s)
        return uniq

    # 回退：拆分
    parts = t.split(delimiter)
    candidates: List[str] = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if "|" in part:
            part = part.split("|", 1)[0].strip()
        if part:
            candidates.append(part)

    # 去重保序
    seen = set()
    out = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def mol_from_smiles_keep_h(smiles: str) -> Optional[Chem.Mol]:
    """
    保留显式 H：禁用 SANITIZE_ADJUSTHS
    """
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None
    try:
        ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
        Chem.SanitizeMol(m, sanitizeOps=ops)
    except Exception:
        return None
    return m


def filter_valid_smiles(candidates: List[str]) -> List[str]:
    valid: List[str] = []
    for s in candidates:
        m = mol_from_smiles_keep_h(s)
        if m is not None:
            valid.append(s)
    return valid


def node_id_to_filename(nid: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", nid.strip())
    if not safe:
        safe = "unknown"
    return f"node_{safe}.png"


def make_mol_from_smiles_list(smiles_list: List[str]) -> Optional[Chem.Mol]:
    if not smiles_list:
        return None

    mols = []
    for s in smiles_list:
        m = mol_from_smiles_keep_h(s)
        if m is None:
            return None
        mols.append(m)

    m = mols[0]
    for m2 in mols[1:]:
        m = Chem.CombineMols(m, m2)

    return m


def prepare_mol_for_drawing(mol: Chem.Mol, explicit_h: bool) -> Chem.Mol:
    """
    兼容新旧 RDKit：
    - 新版：支持 removeHs 参数
    - 旧版：PrepareMolForDrawing 会自动移除 H
      => explicit_h=True 时直接跳过 PrepareMolForDrawing
    """
    m = Chem.Mol(mol)
    if not explicit_h:
        m = Chem.RemoveHs(m)

    rdDepictor.Compute2DCoords(m)

    try:
        # 新版支持 removeHs
        m = rdMolDraw2D.PrepareMolForDrawing(m, removeHs=not explicit_h)
    except Exception:
        # 旧版不支持 removeHs
        if explicit_h:
            return m
        m = rdMolDraw2D.PrepareMolForDrawing(m)

    return m


def make_draw_options(explicit_h: bool) -> rdMolDraw2D.MolDrawOptions:
    opts = rdMolDraw2D.MolDrawOptions()
    if explicit_h:
        # 旧版可能没有 omitTerminalHs
        try:
            opts.omitTerminalHs = False
        except Exception:
            pass
    return opts


def render_mol_to_png(mol: Chem.Mol, size: Tuple[int, int], out_path: Path, explicit_h: bool):
    w, h = size
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    drawer.SetDrawOptions(make_draw_options(explicit_h))
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    out_path.write_bytes(drawer.GetDrawingText())


def render_mol_to_pil(mol: Chem.Mol, size: Tuple[int, int], explicit_h: bool):
    from PIL import Image
    w, h = size
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    drawer.SetDrawOptions(make_draw_options(explicit_h))
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    img = Image.open(BytesIO(png))
    return img.convert("RGB")


def count_h_atoms(mol: Chem.Mol) -> int:
    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")


def draw_per_node_images(
    nodes: List[Tuple[str, Dict[str, str]]],
    key_id: str,
    out_nodes_dir: Path,
    delimiter: str,
    subimg_size: Tuple[int, int],
    explicit_h: bool,
    debug_node: Optional[str],
):
    out_nodes_dir.mkdir(parents=True, exist_ok=True)
    prepared_mols: List[Optional[Chem.Mol]] = []
    legends: List[str] = []

    for nid, data_map in nodes:
        text = data_map.get(key_id, "")

        if debug_node and nid == debug_node:
            print(f"[debug] node {nid} raw text:\n{text}\n")

        candidates = extract_smiles_candidates(text, delimiter=delimiter)

        if debug_node and nid == debug_node:
            print(f"[debug] node {nid} extracted candidates: {candidates}\n")

        smiles_list = filter_valid_smiles(candidates)

        if debug_node and nid == debug_node:
            print(f"[debug] node {nid} valid smiles: {smiles_list}\n")

        mol_raw = make_mol_from_smiles_list(smiles_list)

        out_file = out_nodes_dir / node_id_to_filename(nid)
        if mol_raw is None:
            print(f"[warn] node {nid}: SMILES 解析失败或为空 -> 空白图")
            img = Draw.MolToImage(Chem.MolFromSmiles(""), size=subimg_size)
            img.save(out_file)
            prepared_mols.append(None)
        else:
            if explicit_h and count_h_atoms(mol_raw) == 0:
                print(f"[warn] node {nid}: 显式H开启但分子里H原子数为0 -> {smiles_list}")

            mol_draw = prepare_mol_for_drawing(mol_raw, explicit_h=explicit_h)
            render_mol_to_png(mol_draw, subimg_size, out_file, explicit_h)
            prepared_mols.append(mol_draw)

        legends.append(str(nid))

    return prepared_mols, legends


def normalize_mols_for_grid(mols: List[Optional[Chem.Mol]]) -> List[Chem.Mol]:
    empty = Chem.MolFromSmiles("")
    out: List[Chem.Mol] = []
    for m in mols:
        if isinstance(m, Chem.rdchem.Mol):
            out.append(m)
        else:
            out.append(empty)
    return out


def draw_all_nodes_image(
    mols: List[Optional[Chem.Mol]],
    legends: List[str],
    out_path: Path,
    subimg_size: Tuple[int, int],
    cols: Optional[int],
    explicit_h: bool,
):
    n = len(mols)
    if n == 0:
        raise ValueError("没有节点可绘制")

    if cols is None:
        cols = max(1, int(math.ceil(math.sqrt(n))))

    mols_for_grid = normalize_mols_for_grid(mols)

    try:
        from PIL import Image, ImageDraw
        use_pil = True
    except Exception:
        use_pil = False

    if use_pil:
        cell_w, cell_h = subimg_size
        legend_h = 22
        rows = int(math.ceil(n / cols))
        grid_w = cols * cell_w
        grid_h = rows * (cell_h + legend_h)

        grid = Image.new("RGB", (grid_w, grid_h), "white")
        draw = ImageDraw.Draw(grid)

        for i, mol in enumerate(mols_for_grid):
            r = i // cols
            c = i % cols
            x = c * cell_w
            y = r * (cell_h + legend_h)

            img = render_mol_to_pil(mol, subimg_size, explicit_h)
            grid.paste(img, (x, y))

            if legends:
                draw.text((x + 2, y + cell_h + 2), str(legends[i]), fill="black")

        grid.save(out_path)
        return

    # fallback（如果 PIL 不可用）
    opts = make_draw_options(explicit_h)
    try:
        img = Draw.MolsToGridImage(
            mols_for_grid,
            molsPerRow=cols,
            subImgSize=subimg_size,
            legends=legends,
            useSVG=False,
            prepareMolsBeforeDrawing=False,
            drawOptions=opts,
        )
    except Exception:
        img = Draw.MolsToGridImage(
            mols_for_grid,
            molsPerRow=cols,
            subImgSize=subimg_size,
            legends=legends,
            useSVG=False,
        )
    img.save(out_path)


def main():
    ap = argparse.ArgumentParser("Draw RDKit images from GraphML node SMILES")
    ap.add_argument("--graphml", required=True, help="Input GraphML file")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--smiles-key", default="label_smiles",
                    help="GraphML node key id (e.g. d13) or attr.name (e.g. label_smiles)")
    ap.add_argument("--fallback-keys", default="",
                    help="Comma-separated fallback keys (id or name), used if primary missing")
    ap.add_argument("--delimiter", default=" + ", help="Delimiter between species, default: ' + '")
    ap.add_argument("--cols", type=int, default=None, help="Columns for big image (default: auto sqrt)")
    ap.add_argument("--size", default="300x300", help="Sub-image size, e.g. 300x300")
    ap.add_argument("--explicit-h", action="store_true",
                    help="Draw explicit hydrogens (keep explicit H in SMILES)")
    ap.add_argument("--quiet-rdkit", action="store_true",
                    help="Suppress RDKit parse warnings")
    ap.add_argument("--debug-node", default=None,
                    help="Print raw text and extracted smiles for a given node id")

    args = ap.parse_args()

    if args.quiet_rdkit:
        RDLogger.DisableLog("rdApp.warning")

    size_m = re.match(r"^\s*(\d+)\s*x\s*(\d+)\s*$", args.size)
    if not size_m:
        raise ValueError("size 格式应为 300x300")
    subimg_size = (int(size_m.group(1)), int(size_m.group(2)))

    graphml_path = Path(args.graphml)
    out_dir = Path(args.out_dir)
    out_nodes_dir = out_dir / "nodes"
    out_all = out_dir / "all_nodes.png"

    nodes, key_id_to_name, key_name_to_id = load_graphml_nodes(graphml_path)

    # resolve primary key
    key_id = resolve_key(args.smiles_key, key_id_to_name, key_name_to_id)

    # 若 primary key 在某些节点缺失，尝试 fallback
    fallback_ids: List[str] = []
    if args.fallback_keys.strip():
        for k in args.fallback_keys.split(","):
            k = k.strip()
            if not k:
                continue
            fallback_ids.append(resolve_key(k, key_id_to_name, key_name_to_id))

    # 若节点缺 key，尝试 fallback
    normalized_nodes: List[Tuple[str, Dict[str, str]]] = []
    for nid, data_map in nodes:
        if key_id in data_map and data_map.get(key_id, "").strip():
            normalized_nodes.append((nid, data_map))
            continue
        used = False
        for fk in fallback_ids:
            if fk in data_map and data_map.get(fk, "").strip():
                new_map = dict(data_map)
                new_map[key_id] = data_map[fk]
                normalized_nodes.append((nid, new_map))
                used = True
                break
        if not used:
            normalized_nodes.append((nid, data_map))

    # 排序：若是数字ID则按数字，否则按字符串
    def sort_key(item):
        nid = item[0]
        return (0, int(nid)) if nid.isdigit() else (1, nid)

    normalized_nodes.sort(key=sort_key)

    # 绘制
    out_dir.mkdir(parents=True, exist_ok=True)

    mols, legends = draw_per_node_images(
        normalized_nodes, key_id, out_nodes_dir, args.delimiter, subimg_size,
        args.explicit_h, args.debug_node
    )
    draw_all_nodes_image(mols, legends, out_all, subimg_size, args.cols, args.explicit_h)

    print(f"[done] 单节点图输出至: {out_nodes_dir}")
    print(f"[done] 总览大图: {out_all}")


if __name__ == "__main__":
    # --- IDE quick config ---
    USE_IDE_CONFIG = True

    GRAPHML_IN = "C12H26_network.graphml"
    OUT_DIR = "out_images"

    SMILES_KEY = "label_smiles"   # 也可以写成 "d13"
    FALLBACK_KEYS = "label,orig_id"  # 可选，比如 "d1,d0"

    DELIMITER = " + "
    COLS = None
    SIZE = "300x300"
    EXPLICIT_H = True  # 显式氢（仅保留 SMILES 里已有的）

    if USE_IDE_CONFIG:
        argv = [
            sys.argv[0],
            "--graphml", GRAPHML_IN,
            "--out-dir", OUT_DIR,
            "--smiles-key", SMILES_KEY,
            "--fallback-keys", FALLBACK_KEYS,
            "--delimiter", DELIMITER,
            "--size", SIZE,
        ]
        if COLS is not None:
            argv += ["--cols", str(COLS)]
        if EXPLICIT_H:
            argv += ["--explicit-h"]
        sys.argv = argv
        main()
    else:
        main()