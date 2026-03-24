#!/usr/bin/env python3
VERSION = "104.2"

"""
=============================================================================
PROTEIN MASTER SCREENER (Version 104.2 - INLINE INTERPRETATION)
=============================================================================
Author: Gemini CLI Agent & Kim Junsu, Extended by Claude
Description: 
    Analyzes multi-chain complexes as single entities.
    Each analysis tool outputs to its own folder with individual reports.
    A Summary folder contains the combined report and heatmaps.
=============================================================================
"""

import os
import re
import shutil
import subprocess
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

MIN_SEQ_LENGTH = 5  # Minimum amino acid sequence length to process

D3TO1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
         'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

CONFIG = {
    "INPUT_DIR": "/root/For_phyproper",
    "REF_TAGS": ['test', 'ref', 'WT', '01'],
    "THERMO_ROOT": "/root/ThermoMPNN",
    "GATSOL_ROOT": "/root/GATSol/Predict",
    "MHC_ROOT": "/root/tools",
    "MHC_PY": "/root/mhc1_env/bin/python3",
    "AGGRESCAN_PY": "/root/aggrescan/aggrescan.py",
    "AGGNET_ROOT": "/root/AggNet",
    "NETSOLP_ROOT": "/root/NetSolP/PredictionServer",
    "CONDA": "/root/miniconda3/condabin/conda",
}

# ─── ACCEPTANCE CRITERIA & INTERPRETATION ────────────────────────────────────
# Each metric: (direction, threshold, unit, description)
#   direction: "lower" = lower is better, "higher" = higher is better
#   threshold: acceptance cutoff value
METRIC_CRITERIA = {
    "Stability_ddG": {
        "direction": "lower", "threshold": 0.0, "unit": "kcal/mol",
        "criteria_text": "구조 안정성 변화량(kcal/mol). 음수=안정화 | Criteria: ddG < 0",
        "interpret": lambda v: "PASS - Stabilizing" if v < -0.5 else
                     "PASS - Mildly stabilizing" if v < 0 else
                     "BORDERLINE - Neutral" if v < 0.5 else
                     "FAIL - Destabilizing"
    },
    "Solubility": {
        "direction": "higher", "threshold": 0.5, "unit": "probability",
        "criteria_text": "구조 기반 가용성 예측 확률(0~1). 높을수록 가용성 좋음 | Criteria: > 0.5",
        "interpret": lambda v: "PASS - Highly soluble" if v >= 0.7 else
                     "PASS - Soluble" if v >= 0.5 else
                     "BORDERLINE - Moderate" if v >= 0.4 else
                     "FAIL - Low solubility"
    },
    "Aggrescan_Na4vSS": {
        "direction": "lower", "threshold": 0.0, "unit": "a.u.",
        "criteria_text": "서열 응집 경향 총합 정규화값. 음수=비응집성 | Criteria: < 0",
        "interpret": lambda v: "PASS - Very low aggregation" if v < -10 else
                     "PASS - Low aggregation" if v < 0 else
                     "CAUTION - Mildly aggregation prone" if v < 20 else
                     "FAIL - Highly aggregation prone"
    },
    "Aggrescan_nHS": {
        "direction": "lower", "threshold": 5, "unit": "count",
        "criteria_text": "연속 응집 잔기 구간(핫스팟) 수. 적을수록 좋음 | Criteria: < 5",
        "interpret": lambda v: "PASS - Very few hot spots" if v <= 2 else
                     "PASS - Acceptable" if v <= 5 else
                     "CAUTION - Multiple hot spots" if v <= 15 else
                     "FAIL - Excessive hot spots"
    },
    "Aggrescan_AAT": {
        "direction": "lower", "threshold": 0.0, "unit": "a.u.",
        "criteria_text": "전체 잔기 평균 응집 경향값. 음수=비응집성 | Criteria: < 0",
        "interpret": lambda v: "PASS - Favorable" if v < -0.02 else
                     "PASS - Acceptable" if v < 0 else
                     "CAUTION - Slightly positive" if v < 0.02 else
                     "FAIL - Unfavorable"
    },
    "AggNet_APR_count": {
        "direction": "lower", "threshold": 10, "unit": "residues",
        "criteria_text": "DL 예측 응집성 잔기(APR) 개수. 적을수록 좋음 | Criteria: < 10",
        "interpret": lambda v: "PASS - Minimal APR" if v <= 5 else
                     "PASS - Acceptable" if v <= 10 else
                     "CAUTION - Moderate APR" if v <= 50 else
                     "FAIL - High APR content"
    },
    "AggNet_APR_frac": {
        "direction": "lower", "threshold": 0.05, "unit": "fraction",
        "criteria_text": "전체 서열 중 APR 잔기 비율(0~1). 낮을수록 좋음 | Criteria: < 0.05",
        "interpret": lambda v: "PASS - Very low" if v < 0.03 else
                     "PASS - Acceptable" if v < 0.05 else
                     "CAUTION - Moderate" if v < 0.15 else
                     "FAIL - High APR fraction"
    },
    "AggNet_max_score": {
        "direction": "lower", "threshold": 0.5, "unit": "score",
        "criteria_text": "가장 응집 위험한 잔기의 스코어. 낮을수록 좋음 | Criteria: < 0.5",
        "interpret": lambda v: "PASS - Low peak" if v < 0.3 else
                     "PASS - Acceptable" if v < 0.5 else
                     "CAUTION - Moderate peak" if v < 0.7 else
                     "FAIL - High aggregation peak"
    },
    "NetSolP_Sol": {
        "direction": "higher", "threshold": 0.5, "unit": "probability",
        "criteria_text": "E.coli 발현 시 가용성 확률(0~1). 높을수록 좋음 | Criteria: > 0.5",
        "interpret": lambda v: "PASS - Highly soluble" if v >= 0.7 else
                     "PASS - Soluble" if v >= 0.5 else
                     "CAUTION - Moderate" if v >= 0.3 else
                     "FAIL - Likely insoluble"
    },
    "NetSolP_Usab": {
        "direction": "higher", "threshold": 0.5, "unit": "probability",
        "criteria_text": "E.coli 발현 시 사용 가능 확률(가용성+폴딩). 높을수록 좋음 | Criteria: > 0.5",
        "interpret": lambda v: "PASS - Highly usable" if v >= 0.7 else
                     "PASS - Usable" if v >= 0.5 else
                     "CAUTION - Moderate" if v >= 0.3 else
                     "FAIL - Low usability"
    },
}

# MHC criteria: values are AVERAGE %Rank across ALL sliding-window peptides.
# Most peptides do NOT bind MHC, so average %Rank of 30-50 is typical for antibodies.
# Lower average = more peptides bind MHC = higher immunogenicity risk.
MHC_CRITERIA = {
    "M1": {
        "direction": "higher", "threshold": 50.0, "unit": "avg %Rank",
        "criteria_text": "MHC-I 전체 펩타이드 평균 %Rank. 높을수록 면역원성 낮음 | Criteria: avg > 40",
        "interpret": lambda v: "PASS - Low immunogenicity" if v >= 45 else
                     "PASS - Acceptable" if v >= 40 else
                     "CAUTION - Moderate immunogenicity" if v >= 30 else
                     "FAIL - High immunogenicity risk"
    },
    "M2": {
        "direction": "higher", "threshold": 50.0, "unit": "avg %Rank",
        "criteria_text": "MHC-II 전체 펩타이드 평균 %Rank. 높을수록 면역원성 낮음 | Criteria: avg > 40",
        "interpret": lambda v: "PASS - Low immunogenicity" if v >= 50 else
                     "PASS - Acceptable" if v >= 40 else
                     "CAUTION - Moderate immunogenicity" if v >= 30 else
                     "FAIL - High immunogenicity risk"
    },
}


def get_criteria_for_metric(col_name):
    """Return criteria dict for a given column name."""
    if col_name in METRIC_CRITERIA:
        return METRIC_CRITERIA[col_name]
    if col_name.startswith("M1_"):
        return MHC_CRITERIA["M1"]
    if col_name.startswith("M2_"):
        return MHC_CRITERIA["M2"]
    return None


def add_interpretation(df, metric_columns=None):
    """Merge interpretation into each metric value as 'value (INTERPRETATION)'.
    Also adds a Criteria header row as the first data row for reference.
    Returns a new DataFrame with string-typed metric columns.
    """
    if metric_columns is None:
        metric_columns = [c for c in df.columns if c not in ('ID', 'Chains', 'Chain_detail', 'Total_length')]

    df_out = df.copy()

    # Build criteria row (first row showing acceptance criteria for each metric)
    criteria_row = {}
    for c in df_out.columns:
        if c in metric_columns:
            crit = get_criteria_for_metric(c)
            criteria_row[c] = f"[Criteria: {crit['criteria_text']}]" if crit else ""
        elif c == 'ID':
            criteria_row[c] = '--- CRITERIA ---'
        elif c == 'Chains':
            criteria_row[c] = '분석된 체인 ID'
        elif c == 'Chain_detail':
            criteria_row[c] = '체인별 잔기 수'
        elif c == 'Total_length':
            criteria_row[c] = '총 아미노산 수'
        else:
            criteria_row[c] = ''

    # Merge interpretation into each cell: "value (INTERPRETATION)"
    for c in metric_columns:
        if c not in df_out.columns:
            continue
        crit = get_criteria_for_metric(c)
        if crit:
            def _safe_int_check(v):
                if pd.isna(v): return False
                return isinstance(v, int) or (isinstance(v, float) and v == int(v))
            is_int_col = all(_safe_int_check(v) for v in df_out[c])
            def _fmt(v, is_int, crit):
                if pd.isna(v): return "N/A (Tool Failed)"
                if is_int:
                    return f"{int(v)} ({crit['interpret'](v)})"
                return f"{v:.4f} ({crit['interpret'](v)})"
            df_out[c] = [_fmt(v, is_int_col, crit) for v in df_out[c]]

    # Insert criteria row at the top
    criteria_df = pd.DataFrame([criteria_row], columns=df_out.columns)
    df_out = pd.concat([criteria_df, df_out], ignore_index=True)

    return df_out



class AssemblySelect(Select):
    def accept_residue(self, residue): return 1 if residue.id[0] == " " else 0

_log_file = None  # Global log file path, set in main()

def _log(msg):
    """Write a message to the reproducibility log file."""
    if _log_file and _log_file.exists():
        with open(_log_file, 'a') as f:
            f.write(msg + "\n")

def run_cmd(cmd, cwd=None, timeout=600):
    """Run a command and return stdout. Returns None on failure."""
    cmd_str = " ".join(str(c) for c in cmd)
    _log(f"[EXEC] {cmd_str}")
    try:
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            print(f"    [WARN] Command failed (rc={result.returncode}): {cmd_str}")
            _log(f"[FAIL] rc={result.returncode}")
            if result.stderr:
                print(f"    stderr: {result.stderr[:300]}")
            return None
        return result.stdout
    except subprocess.TimeoutExpired:
        print(f"    [WARN] Command timed out after {timeout}s")
        return None
    except Exception as e:
        print(f"    [WARN] Command error: {e}")
        return None

def get_standardized_assembly(input_p, output_p):
    """Parse structure, renumber residues globally, return (full_seq, chain_info_list).
    chain_info_list: [{"chain_id": "A", "length": 107, "seq": "ETTVT..."}, ...]
    """
    parser = MMCIFParser(QUIET=True) if input_p.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    struct = parser.get_structure('x', str(input_p))
    full_seq_list = []
    chain_info_list = []
    model = struct[0]  # Process only the first model (C5 fix)
    if True:  # preserving indentation block
        # Collect original chain IDs and sequences before renumbering
        orig_chains = []
        for chain in model:
            chain_seq = []
            for residue in chain.get_residues():
                if residue.id[0] == " ":
                    aa = D3TO1.get(residue.get_resname(), '')
                    if aa: chain_seq.append(aa)
            if chain_seq:
                orig_chains.append({"chain_id": chain.id, "seq": "".join(chain_seq)})
        # Renumber residues globally
        temp_id = 90000
        for chain in model:
            for residue in list(chain.get_residues()):
                if residue.id[0] == " ":
                    temp_id += 1
                    residue.id = (' ', temp_id, ' ')
        global_idx = 0
        for chain in model:
            chain_seq = []
            for residue in list(chain.get_residues()):
                if residue.id[0] == " ":
                    global_idx += 1
                    residue.id = (' ', global_idx, ' ')
                    aa = D3TO1.get(residue.get_resname(), '')
                    if aa: chain_seq.append(aa)
            if chain_seq: full_seq_list.append("".join(chain_seq))
        # Build chain_info from original data
        for ci in orig_chains:
            chain_info_list.append({
                "chain_id": ci["chain_id"],
                "length": len(ci["seq"]),
                "seq": ci["seq"],
            })
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(output_p), AssemblySelect())
    return "".join(full_seq_list), chain_info_list


# ─── 1. ThermoMPNN ──────────────────────────────────────────────────────────
def run_thermompnn(work_pdb, stem, tool_dir):
    """Run ThermoMPNN stability prediction. Returns dict and raw output."""
    out_s = run_cmd([CONFIG["CONDA"], "run", "-n", "thermoMPNN", "--cwd", CONFIG["THERMO_ROOT"],
                     "python3", "minimal_ddg.py", "--pdb_path", str(work_pdb)])
    stab = float('nan')
    if out_s is not None:
        m = re.search(r"Predicted .{1,5}G \(kcal/mol\):\s*([-0-9.]+)", out_s)
        if m:
            val = float(m.group(1))
            if abs(val) <= 100:  # F8: sanity check
                stab = val
            else:
                print(f"    [WARN] Implausible ddG value: {val}")
        # Save raw output
        with open(tool_dir / f"{stem}_raw_output.txt", 'w') as f:
            f.write(out_s)
    return {"Stability_ddG": stab}, out_s or ""


# ─── 2. GATSol ──────────────────────────────────────────────────────────────
def run_gatsol(work_pdb, full_seq, stem, tool_dir):
    """Run GATSol solubility prediction. Returns dict."""
    g_root = Path(CONFIG["GATSOL_ROOT"])
    for s in ["pdb", "fasta", "cm", "pkl", "output"]:
        p = g_root / "NEED_to_PREPARE" / s
        if p.exists():
            for f in p.glob("*"):
                if f.is_file():
                    f.unlink()
    shutil.rmtree(g_root / "output", ignore_errors=True)
    os.makedirs(g_root / "output", exist_ok=True)
    shutil.copy(work_pdb, g_root / "NEED_to_PREPARE/pdb" / f"{stem}.pdb")
    with open(g_root / "NEED_to_PREPARE/list.csv", "w") as f:
        f.write(f"id,sequence\n{stem},{full_seq}\n")
    with open(g_root / "NEED_to_PREPARE/fasta" / f"{stem}.fasta", "w") as f:
        f.write(f">{stem}\n{full_seq}\n")
    run_cmd([CONFIG["CONDA"], "run", "-n", "GATSol", "--cwd", str(g_root), "bash", "./tools/Predict.sh"])
    sol = float('nan')
    sol_csv = g_root / "output/prediction_results.csv"
    if sol_csv.exists():
        try:
            df_sol = pd.read_csv(sol_csv, dtype={'id': str})
            match = df_sol[df_sol['id'] == stem]
            sol = float(match.iloc[0, 1]) if not match.empty else float(df_sol.iloc[0, 1])
            shutil.copy(sol_csv, tool_dir / f"{stem}_prediction_results.csv")
        except Exception as e:
            print(f"    [WARN] GATSol result parsing failed: {e}")
    return {"Solubility": sol}


# ─── 3. MHC Pipeline ────────────────────────────────────────────────────────
def run_mhc(full_seq, stem, tool_dir):
    """Run MHC I/II immunogenicity prediction. Returns dict."""
    mhc_out = Path(CONFIG["MHC_ROOT"]) / "output"
    shutil.rmtree(mhc_out, ignore_errors=True)
    os.makedirs(mhc_out, exist_ok=True)
    fa_m = Path(CONFIG["MHC_ROOT"]) / f"tmp_{stem}.fasta"
    with open(fa_m, "w") as f:
        f.write(f">{stem}\n{full_seq}\n")
    run_cmd([CONFIG["MHC_PY"], "mhc_pipeline.py", "--input", str(fa_m)], cwd=CONFIG["MHC_ROOT"])
    if fa_m.exists():
        fa_m.unlink()  # cleanup temp fasta
    mhc_row = {}
    for csv_file in (Path(CONFIG["MHC_ROOT"]) / "output").rglob("combined_MHC*.csv"):
        dm = pd.read_csv(csv_file)
        pfx = "M1_" if "MHC1" in csv_file.name else "M2_"
        for col in dm.columns[1:]:
            mhc_row[f"{pfx}{col}"] = dm.iloc[0][col]
        # Copy raw CSV to tool folder
        shutil.copy(csv_file, tool_dir / f"{stem}_{csv_file.name}")
    return mhc_row


# ─── 4. AGGRESCAN ───────────────────────────────────────────────────────────
def run_aggrescan(seq, stem, tool_dir):
    """Run AGGRESCAN sequence-based aggregation prediction."""
    try:
        run_cmd(["python3", CONFIG["AGGRESCAN_PY"], "-s", seq, "-n", stem,
                 "-o", str(tool_dir), "--no-plot"])
        summary = tool_dir / "aggrescan_summary.tsv"
        if summary.exists():
            df = pd.read_csv(summary, sep='\t')
            row = df[df['id'] == stem]
            if not row.empty:
                return {
                    "Aggrescan_Na4vSS": float(row.iloc[0]['Na4vSS']),
                    "Aggrescan_nHS": int(row.iloc[0]['nHS']),
                    "Aggrescan_AAT": float(row.iloc[0]['AAT']),
                }
    except Exception as e:
        print(f"    [WARN] AGGRESCAN failed: {e}")
    return {"Aggrescan_Na4vSS": float("nan"), "Aggrescan_nHS": float("nan"), "Aggrescan_AAT": float("nan")}


# ─── 5. AggNet ──────────────────────────────────────────────────────────────
def run_aggnet(seq, stem, tool_dir):
    """Run AggNet DL-based aggregation-prone region prediction."""
    aggnet_out = tool_dir / f"{stem}_per_residue.csv"
    try:
        run_cmd([CONFIG["CONDA"], "run", "-n", "AggNet", "--cwd", CONFIG["AGGNET_ROOT"],
                 "python", "./script/predict_APR.py",
                 "--sequence", seq, "--checkpoint", "./checkpoint/APNet.ckpt",
                 "--output", str(aggnet_out)])
        if aggnet_out.exists():
            df = pd.read_csv(aggnet_out)
            n_total = len(df)
            n_apr = int(df['APR'].sum())
            return {
                "AggNet_APR_count": n_apr,
                "AggNet_APR_frac": round(n_apr / n_total, 4) if n_total > 0 else 0.0,
                "AggNet_max_score": round(float(df['scores'].max()), 4),
            }
    except Exception as e:
        print(f"    [WARN] AggNet failed: {e}")
    return {"AggNet_APR_count": float("nan"), "AggNet_APR_frac": float("nan"), "AggNet_max_score": float("nan")}


# ─── 6. NetSolP ─────────────────────────────────────────────────────────────
def run_netsolp(seq, stem, tool_dir):
    """Run NetSolP E. coli solubility/usability prediction."""
    fasta_p = tool_dir / f"{stem}.fasta"
    csv_p = tool_dir / f"{stem}_prediction.csv"
    try:
        with open(fasta_p, 'w') as f:
            f.write(f">{stem}\n{seq}\n")
        run_cmd([CONFIG["CONDA"], "run", "-n", "netsolp", "--cwd", CONFIG["NETSOLP_ROOT"],
                 "python", "predict.py",
                 "--FASTA_PATH", str(fasta_p), "--OUTPUT_PATH", str(csv_p),
                 "--MODEL_TYPE", "ESM12", "--PREDICTION_TYPE", "SU"])
        if csv_p.exists():
            df = pd.read_csv(csv_p)
            if not df.empty:
                return {
                    "NetSolP_Sol": round(float(df.iloc[0]['predicted_solubility']), 4),
                    "NetSolP_Usab": round(float(df.iloc[0]['predicted_usability']), 4),
                }
    except Exception as e:
        print(f"    [WARN] NetSolP failed: {e}")
    return {"NetSolP_Sol": float("nan"), "NetSolP_Usab": float("nan")}


# ─── HEATMAP GENERATION ─────────────────────────────────────────────────────
def generate_heatmaps(df_f, ref_id, summary_dir):
    """Generate comparison and safety-zone heatmaps in Summary/Visualizations/."""
    # Drop non-numeric info columns before heatmap
    drop_cols = [c for c in df_f.columns if c in ('Chains', 'Chain_detail', 'Total_length')]
    df_p = df_f.drop(columns=drop_cols, errors='ignore').set_index('ID')
    # Fill NaN with 0 for heatmap rendering (NaN would break seaborn)
    df_p = df_p.fillna(0.0)
    viz_dir = summary_dir / "Visualizations"
    viz_dir.mkdir(exist_ok=True)

    # Derive direction from METRIC_CRITERIA (single source of truth)
    higher_better = {k: v["direction"] == "higher" for k, v in METRIC_CRITERIA.items()}
    for c in df_p.columns:
        if c.startswith("M1_"):
            higher_better[c] = MHC_CRITERIA["M1"]["direction"] == "higher"
        elif c.startswith("M2_"):
            higher_better[c] = MHC_CRITERIA["M2"]["direction"] == "higher"

    # 1. Comparison Heatmap (vs Reference)
    if ref_id and ref_id in df_p.index:
        ref = df_p.loc[ref_id]
        diff_wt = pd.DataFrame(index=df_p.index, columns=df_p.columns, dtype=float)
        for c in df_p.columns:
            if higher_better.get(c, False):
                diff_wt[c] = df_p[c] - ref[c]
            else:
                diff_wt[c] = ref[c] - df_p[c]
        fig_w = max(28, len(diff_wt.columns) * 1.8)
        plt.figure(figsize=(fig_w, max(6, len(diff_wt) * 1.2)))
        sns.heatmap(diff_wt.astype(float), annot=True, cmap='RdYlBu', center=0, fmt='.2f',
                    linewidths=0.5, linecolor='gray')
        plt.title(f'IMPROVEMENT vs {ref_id} (Blue: Better)', fontsize=20)
        plt.tight_layout()
        plt.savefig(viz_dir / "1_Comparison_vs_WT.png", dpi=300)
        plt.close()

    # 2. Safety Zone Heatmap - derive thresholds from METRIC_CRITERIA
    safety_thresholds = {k: (v["direction"], float(v["threshold"])) for k, v in METRIC_CRITERIA.items()}
    for c in df_p.columns:
        if c.startswith("M1_"):
            safety_thresholds[c] = (MHC_CRITERIA["M1"]["direction"], float(MHC_CRITERIA["M1"]["threshold"]))
        elif c.startswith("M2_"):
            safety_thresholds[c] = (MHC_CRITERIA["M2"]["direction"], float(MHC_CRITERIA["M2"]["threshold"]))

    diff_s = pd.DataFrame(index=df_p.index, columns=df_p.columns, dtype=float)
    for c in df_p.columns:
        direction, threshold = safety_thresholds.get(c, ('lower', 0.0))
        if direction == 'higher':
            diff_s[c] = df_p[c] - threshold
        else:
            diff_s[c] = threshold - df_p[c]
    fig_w = max(28, len(diff_s.columns) * 1.8)
    plt.figure(figsize=(fig_w, max(6, len(diff_s) * 1.2)))
    # Ensure shapes match (C2 fix)
    assert diff_s.shape == df_p.shape, f"Shape mismatch: diff_s={diff_s.shape}, df_p={df_p.shape}"
    df_p_numeric = df_p.apply(pd.to_numeric, errors='coerce')  # F9: ensure numeric annotations
    sns.heatmap(diff_s.astype(float), annot=df_p_numeric, cmap='RdYlBu', center=0, fmt='.2f',
                linewidths=0.5, linecolor='gray')
    plt.title('SAFETY ZONE ANALYSIS (Blue: Safe Zone)', fontsize=20)
    plt.tight_layout()
    plt.savefig(viz_dir / "2_Safety_Threshold_Analysis.png", dpi=300)
    plt.close()


# ─── MAIN ────────────────────────────────────────────────────────────────────
def main():
    in_path = Path(CONFIG["INPUT_DIR"])
    if not in_path.exists():
        print(f"[ERROR] Input directory not found: {in_path}")
        return
    # Date-based output directory with auto-increment
    from datetime import datetime
    result_base = in_path / "result"
    today = datetime.now().strftime("%Y-%m-%d")
    out_root = result_base / today
    if out_root.exists():
        run_num = 1
        while (result_base / f"{today}_{run_num}").exists():
            run_num += 1
        out_root = result_base / f"{today}_{run_num}"
    out_root.mkdir(parents=True, exist_ok=True)

    # Create per-tool directories
    dirs = {
        "thermo":    out_root / "01_ThermoMPNN",
        "gatsol":    out_root / "02_GATSol",
        "mhc":       out_root / "03_MHC",
        "aggrescan": out_root / "04_AGGRESCAN",
        "aggnet":    out_root / "05_AggNet",
        "netsolp":   out_root / "06_NetSolP",
        "summary":   out_root / "Summary",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # ── Reproducibility Log ──
    global _log_file
    _log_file = out_root / "pipeline_run.log"
    import platform
    with open(_log_file, 'w') as lf:
        lf.write("=" * 70 + "\n")
        lf.write("  Protein Master Screener - Reproducibility Log\n")
        lf.write("=" * 70 + "\n\n")
        lf.write(f"Date/Time       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        lf.write(f"Version         : {VERSION}\n")
        lf.write(f"Input Dir       : {in_path}\n")
        lf.write(f"Output Dir      : {out_root}\n")
        lf.write(f"Python          : {__import__('sys').version.split()[0]}\n")
        lf.write(f"Platform        : {platform.platform()}\n\n")
        lf.write("── Configuration ──\n")
        for k, v in CONFIG.items():
            lf.write(f"  {k:15s}: {v}\n")
        lf.write("\n── Tool Versions ──\n")
    # ── Tool Versions (system + per-conda-env) ──
    for tool_cmd, label in [
        (["python3", "--version"], "Python3"),
        ([CONFIG["CONDA"], "--version"], "Conda"),
    ]:
        try:
            r = subprocess.run(tool_cmd, capture_output=True, text=True, timeout=10)
            ver = (r.stdout.strip() or r.stderr.strip()).split("\n")[0]
            _log(f"  {label:15s}: {ver}")
        except Exception:
            _log(f"  {label:15s}: (not found)")
    # System Python libraries
    for mod_name, label in [("Bio", "BioPython"), ("pandas", "Pandas")]:
        try:
            mod = __import__(mod_name)
            _log(f"  {label:15s}: {getattr(mod, '__version__', 'unknown')}")
        except ImportError:
            _log(f"  {label:15s}: (not installed)")
    # Per-tool conda environment versions
    _log("")
    _log("  [Per-tool environment versions]")
    env_checks = [
        ("thermoMPNN", "import torch; import numpy; import omegaconf; "
         "print(f'PyTorch={torch.__version__}, NumPy={numpy.__version__}, OmegaConf={omegaconf.__version__}')"),
        ("GATSol", "import torch; import torch_geometric; print(f'PyTorch={torch.__version__}, PyG={torch_geometric.__version__}')"),
        ("AggNet", "import torch; import pytorch_lightning as pl; "
         "print(f'PyTorch={torch.__version__}, PL={pl.__version__}')"),
        ("netsolp", "import onnxruntime; print(f'ONNX-RT={onnxruntime.__version__}')"),
    ]
    for env_name, py_code in env_checks:
        try:
            r = subprocess.run(
                [CONFIG["CONDA"], "run", "-n", env_name, "python3", "-c", py_code],
                capture_output=True, text=True, timeout=30)
            ver = (r.stdout.strip() or r.stderr.strip()).split("\n")[-1] if r.returncode == 0 else "(env error)"
            _log(f"    {env_name:14s}: {ver}")
        except Exception:
            _log(f"    {env_name:14s}: (not available)")
    # MHC environment
    try:
        r = subprocess.run([CONFIG["MHC_PY"], "-c", "import sys; print(f'Python={sys.version.split()[0]}')"],
                           capture_output=True, text=True, timeout=10)
        _log(f"    {'mhc1_env':14s}: {r.stdout.strip()}")
    except Exception:
        _log(f"    {'mhc1_env':14s}: (not available)")

    # Per-tool software versions (git commit or tag)
    _log("")
    _log("  [Tool software versions (git)]")
    tool_repos = [
        ("ThermoMPNN", CONFIG["THERMO_ROOT"]),
        ("GATSol",     str(Path(CONFIG["GATSOL_ROOT"]).parent)),
        ("AggNet",     CONFIG["AGGNET_ROOT"]),
        ("NetSolP",    str(Path(CONFIG["NETSOLP_ROOT"]).parent)),
    ]
    for tool_label, repo_path in tool_repos:
        try:
            # Try git describe --tags first, fall back to commit hash
            tag_r = subprocess.run(["git", "-C", repo_path, "describe", "--tags"],
                                   capture_output=True, text=True, timeout=5)
            if tag_r.returncode == 0 and tag_r.stdout.strip():
                ver_str = tag_r.stdout.strip()
            else:
                hash_r = subprocess.run(["git", "-C", repo_path, "log", "--oneline", "-1",
                                         "--format=%h %ci"], capture_output=True, text=True, timeout=5)
                ver_str = hash_r.stdout.strip() if hash_r.returncode == 0 else "(unknown)"
            origin_r = subprocess.run(["git", "-C", repo_path, "remote", "get-url", "origin"],
                                      capture_output=True, text=True, timeout=5)
            origin = origin_r.stdout.strip() if origin_r.returncode == 0 else ""
            # Check for local modifications
            dirty_r = subprocess.run(["git", "-C", repo_path, "status", "--porcelain"],
                                     capture_output=True, text=True, timeout=5)
            dirty = " [LOCAL MODIFIED]" if (dirty_r.returncode == 0 and dirty_r.stdout.strip()) else ""
            _log(f"    {tool_label:14s}: {ver_str}{dirty}  ({origin})")
        except Exception:
            _log(f"    {tool_label:14s}: (not a git repo)")
    # MixMHCpred / MixMHC2pred (not git repos, check binary/script)
    # MixMHCpred version (from -h output first line)
    try:
        r = subprocess.run([str(Path(CONFIG["MHC_ROOT"]) / "MixMHCpred" / "MixMHCpred"), "-h"],
                           capture_output=True, text=True, timeout=10)
        ver = (r.stdout.strip().split("\n")[0]) if r.stdout.strip() else "(unknown)"
        _log(f"    {'MixMHCpred':14s}: {ver}")
    except Exception:
        _log(f"    {'MixMHCpred':14s}: (not available)")
    # MixMHC2pred version (from NEWS.md)
    try:
        news_path = Path(CONFIG["MHC_ROOT"]) / "MixMHC2pred" / "NEWS.md"
        if news_path.exists():
            import re as _re
            with open(news_path) as _nf:
                for line in _nf:
                    m = _re.search(r"Version\s+([\d.]+)", line)
                    if m:
                        _log(f"    {'MixMHC2pred':14s}: {m.group(1)}")
                        break
                else:
                    _log(f"    {'MixMHC2pred':14s}: (version not found in NEWS.md)")
        else:
            _log(f"    {'MixMHC2pred':14s}: (NEWS.md not found)")
    except Exception:
        _log(f"    {'MixMHC2pred':14s}: (not available)")
    # AGGRESCAN (local implementation based on published algorithm)
    _log(f"    {'AGGRESCAN':14s}: local impl (de Groot 2005 / Conchillo-Sole 2007)")

    _log("")
    _log("── Tool Configurations (applied options per program) ──")
    _log("")
    _log("  [1/6] ThermoMPNN ─ Protein Stability (ddG)")
    _log(f"    Conda env       : thermoMPNN")
    _log(f"    Model weights   : {CONFIG['THERMO_ROOT']}/models/thermoMPNN_default.pt")
    _log(f"    Architecture    : hidden_dims=[64,32], num_final_layers=2, lightattn=True")
    _log(f"    ProteinMPNN     : model_type=vanilla, freeze_weights=True, load_pretrained=True")
    _log(f"    Training cfg    : subtract_mut=True")
    _log(f"    Random seed     : 42 (deterministic=True, benchmark=False)")
    _log(f"    Device          : CUDA if available, else CPU")
    _log(f"    Applied command : conda run -n thermoMPNN python3 minimal_ddg.py --pdb_path <pdb>")
    _log(f"    Validity guard  : |ddG| <= 100 kcal/mol")
    _log("")
    _log("  [2/6] GATSol ─ Structure-based Solubility")
    _log(f"    Conda env       : GATSol")
    _log(f"    Model weights   : {CONFIG['GATSOL_ROOT']}/../check_point/best_model/best_model.pt")
    _log(f"    Architecture    : GAT 2-layer, hidden=1024, heads=16, pooling=global_mean → FC(128) → 1")
    _log(f"    Node features   : 1300 dims (ESM-1b 1280d layer33 + BLOSUM62 20d)")
    _log(f"    ESM model       : esm1b_t33_650M_UR50S")
    _log(f"    Contact map     : CA atoms, distance threshold=10.0 Å, self-loops=True")
    _log(f"    Batch size      : 1 (inference)")
    _log(f"    Applied command : conda run -n GATSol bash ./tools/Predict.sh")
    _log("")
    _log("  [3/6] MHC I/II ─ Immunogenicity (MixMHCpred / MixMHC2pred)")
    _log(f"    Python env      : {CONFIG['MHC_PY']}")
    _log(f"    MHC-I tool      : {CONFIG['MHC_ROOT']}/MixMHCpred/code/main.py")
    _log(f"    MHC-I alleles   : HLA-A01:01, A02:01, A03:01, B07:02, B08:01, B15:01, C07:01, C04:01")
    _log(f"    MHC-I peptide   : k-mer lengths = 9, 10, 11")
    _log(f"    MHC-II tool     : {CONFIG['MHC_ROOT']}/MixMHC2pred/MixMHC2pred_unix")
    _log(f"    MHC-II alleles  : DRB1_01_01, 03_01, 04_01, 07_01, 09_01, 11_01, 13_01, 15_01")
    _log(f"    MHC-II peptide  : k-mer length = 15, flag = --no_context")
    _log(f"    Scoring         : mean %Rank across all sliding-window peptides per allele")
    _log(f"    Default score   : 100.0 (if tool output missing)")
    _log(f"    Applied command : {CONFIG['MHC_PY']} mhc_pipeline.py --input <fasta>")
    _log("")
    _log("  [4/6] AGGRESCAN ─ Sequence Aggregation (classic)")
    _log(f"    Python          : system python3")
    _log(f"    Scale           : a3v (de Groot et al. 2005, 20 amino acid propensity values)")
    _log(f"    Window size     : adaptive (≤75aa→5, 76-175→7, 176-300→9, >300→11)")
    _log(f"    Hot-spot threshold (HST): ~-0.0929 (SwissProt frequency-weighted a3v average)")
    _log(f"    Hot-spot rule   : a4v > HST AND residue ≠ Pro, min region length = 5")
    _log(f"    Output metrics  : Na4vSS (total), AAT (avg), nHS (hot-spot count)")
    _log(f"    Applied command : python3 aggrescan.py -s <seq> -n <id> -o <dir> --no-plot")
    _log("")
    _log("  [5/6] AggNet ─ Aggregation-Prone Regions (Deep Learning)")
    _log(f"    Conda env       : AggNet")
    _log(f"    Checkpoint      : {CONFIG['AGGNET_ROOT']}/checkpoint/APNet.ckpt")
    _log(f"    Seq-only params : t_start=0.46, t_expand=0.37, t_patience=7")
    _log(f"    Struct params   : beta=3.36, delta=0.4, t_start=0.51, t_expand=0.37, t_patience=9")
    _log(f"    Mode used       : sequence-only (no PDB input)")
    _log(f"    Output metrics  : APR count, APR fraction, max aggregation score")
    _log(f"    Applied command : conda run -n AggNet python predict_APR.py --sequence <seq> --checkpoint <ckpt> --output <csv>")
    _log("")
    _log("  [6/6] NetSolP ─ E.coli Solubility/Usability")
    _log(f"    Conda env       : netsolp")
    _log(f"    Model type      : ESM12 (5-fold cross-validation, quantized ONNX)")
    _log(f"    Prediction type : SU (Solubility + Usability)")
    _log(f"    Activation      : Sigmoid on raw predictions")
    _log(f"    Max seq length  : 1022 residues")
    _log(f"    Applied command : conda run -n netsolp python predict.py --FASTA_PATH <fa> --OUTPUT_PATH <csv> --MODEL_TYPE ESM12 --PREDICTION_TYPE SU")
    _log("")
    _log("  [Shared] Structure Preprocessing")
    _log(f"    Parser          : BioPython MMCIFParser / PDBParser (QUIET=True)")
    _log(f"    Model select    : First model only (model index 0)")
    _log(f"    Residue filter  : Standard residues only (HETATM excluded)")
    _log(f"    Renumbering     : Global sequential (1, 2, 3, ...) across all chains")
    _log(f"    Min seq length  : {MIN_SEQ_LENGTH} aa")
    _log("")
    _log("  [Shared] Acceptance Criteria")
    for metric_name, mc in METRIC_CRITERIA.items():
        _log(f"    {metric_name:22s}: {mc['direction']:>6s} is better, threshold={mc['threshold']} {mc['unit']}")
    for mhc_key, mc in MHC_CRITERIA.items():
        _log(f"    {mhc_key:22s}: {mc['direction']:>6s} is better, threshold={mc['threshold']} {mc['unit']}")


    _log("── Input Files ──")

    raw_files = sorted(list(in_path.glob("*.pdb")) + list(in_path.glob("*.cif")))
    stems = [f.stem for f in raw_files]
    # Use word-boundary matching to avoid false positives (e.g. "contested" matching "test")
    def _match_ref_tag(tag, stem):
        return re.search(rf'(?:^|[_\-])({re.escape(tag)})(?:$|[_\-])', stem, re.IGNORECASE) or stem.lower() == tag.lower()
    ref_id = next((s for t in CONFIG["REF_TAGS"] for s in stems if _match_ref_tag(t, s)), None)
    for rf in raw_files:
        _log(f"  {rf.name}")
    _log(f"  Reference: {ref_id}")
    _log("")
    _log("── Execution Log ──")
    print(f"Assembly Engine v{VERSION} started. Reference: {ref_id}")
    print(f"Output: {out_root}")

    # Accumulators for per-tool reports
    all_thermo, all_gatsol, all_mhc = [], [], []
    all_aggrescan, all_aggnet, all_netsolp = [], [], []
    final_data = []

    for raw_f in raw_files:
        stem = raw_f.stem
        print(f"\n{'='*60}\n  Analyzing Assembly: {stem}\n{'='*60}")
        _log(f"\n[SAMPLE] {stem}")
        work_pdb = out_root / f"{stem}.pdb"
        full_seq, chain_info = get_standardized_assembly(raw_f, work_pdb)
        if len(full_seq) < MIN_SEQ_LENGTH:
            print(f"    [SKIP] Sequence too short ({len(full_seq)} aa, min={MIN_SEQ_LENGTH})")
            continue

        # Chain summary
        chains_str = "+".join([f"{c['chain_id']}({c['length']})" for c in chain_info])
        chain_ids = "+".join([c['chain_id'] for c in chain_info])
        total_len = sum(c['length'] for c in chain_info)
        print(f"    Chains: {chains_str}  Total: {total_len} aa")

        # 1. ThermoMPNN
        print("    [1/6] ThermoMPNN (Stability)...")
        thermo_row, _ = run_thermompnn(work_pdb, stem, dirs["thermo"])
        all_thermo.append({"ID": stem, **thermo_row})

        # 2. GATSol
        print("    [2/6] GATSol (Solubility)...")
        gatsol_row = run_gatsol(work_pdb, full_seq, stem, dirs["gatsol"])
        all_gatsol.append({"ID": stem, **gatsol_row})

        # 3. MHC
        print("    [3/6] MHC Pipeline (Immunogenicity)...")
        mhc_row = run_mhc(full_seq, stem, dirs["mhc"])
        all_mhc.append({"ID": stem, **mhc_row})

        # 4. AGGRESCAN
        print("    [4/6] AGGRESCAN (Aggregation - classic)...")
        aggrescan_row = run_aggrescan(full_seq, stem, dirs["aggrescan"])
        all_aggrescan.append({"ID": stem, **aggrescan_row})

        # 5. AggNet
        print("    [5/6] AggNet (Aggregation - DL SOTA)...")
        aggnet_row = run_aggnet(full_seq, stem, dirs["aggnet"])
        all_aggnet.append({"ID": stem, **aggnet_row})

        # 6. NetSolP
        print("    [6/6] NetSolP (E.coli Solubility)...")
        netsolp_row = run_netsolp(full_seq, stem, dirs["netsolp"])
        all_netsolp.append({"ID": stem, **netsolp_row})

        # Combined row
        combined = {
            "ID": stem,
            "Chains": chain_ids,
            "Chain_detail": chains_str,
            "Total_length": total_len,
            **thermo_row,
            **gatsol_row,
            **mhc_row,
            **aggrescan_row,
            **aggnet_row,
            **netsolp_row,
        }
        final_data.append(combined)

        def _fmt_val(v, fmt=".2f"):
            return f"{v:{fmt}}" if not pd.isna(v) else "N/A"
        print(f"    Done: {stem} | ddG={_fmt_val(thermo_row['Stability_ddG'])} Sol={_fmt_val(gatsol_row['Solubility'])} "
              f"nHS={_fmt_val(aggrescan_row['Aggrescan_nHS'], '.0f')} APR={_fmt_val(aggnet_row['AggNet_APR_count'], '.0f')} "
              f"NetSolP={_fmt_val(netsolp_row['NetSolP_Sol'], '.3f')}")

    if not final_data:
        print("\n[ERROR] No valid samples processed.")
        return

    # ─── Save per-tool reports (with interpretation) ───────────────────
    print(f"\n{'='*60}\n  Saving per-tool reports...\n{'='*60}")

    tool_reports = [
        ("thermo",    all_thermo,    "ThermoMPNN_report.tsv",  "01_ThermoMPNN"),
        ("gatsol",    all_gatsol,     "GATSol_report.tsv",      "02_GATSol"),
        ("mhc",       all_mhc,       "MHC_report.tsv",         "03_MHC"),
        ("aggrescan", all_aggrescan, "AGGRESCAN_report.tsv",    "04_AGGRESCAN"),
        ("aggnet",    all_aggnet,    "AggNet_report.tsv",       "05_AggNet"),
        ("netsolp",   all_netsolp,   "NetSolP_report.tsv",     "06_NetSolP"),
    ]
    for key, data_list, fname, label in tool_reports:
        df_tool = pd.DataFrame(data_list)
        df_tool_interp = add_interpretation(df_tool.copy())
        df_tool_interp.to_csv(dirs[key] / fname, sep='\t', index=False)
        print(f"    {label}/{fname}")

    # ─── Save Summary report ────────────────────────────────────────────
    print(f"\n{'='*60}\n  Generating Summary report & heatmaps...\n{'='*60}")

    df_f = pd.DataFrame(final_data)
    # Save raw data (for heatmaps)
    df_f.to_csv(dirs["summary"] / "Final_Assembly_Report_raw.tsv", sep='\t', index=False)
    # Save with interpretation
    df_f_interp = add_interpretation(df_f.copy())
    df_f_interp.to_excel(dirs["summary"] / "Final_Assembly_Report.xlsx", index=False)
    df_f_interp.to_csv(dirs["summary"] / "Final_Assembly_Report.tsv", sep='\t', index=False)
    print(f"    Summary/Final_Assembly_Report.xlsx (with interpretation)")
    print(f"    Summary/Final_Assembly_Report.tsv (with interpretation)")
    print(f"    Summary/Final_Assembly_Report_raw.tsv (numeric only)")

    generate_heatmaps(df_f, ref_id, dirs["summary"])
    print(f"    Summary/Visualizations/ generated")

    # ── Finalize Log ──
    end_time = datetime.now()
    _log("")
    _log("── Results Summary ──")
    _log(f"  Output Dir    : {out_root}")
    _log(f"  Summary Report: {dirs['summary'] / 'Final_Assembly_Report.xlsx'}")
    _log(f"  Raw Data      : {dirs['summary'] / 'Final_Assembly_Report_raw.tsv'}")
    _log(f"  Heatmaps      : {dirs['summary'] / 'Visualizations/'}")
    _log("")
    _log("── Reproduce This Run ──")
    _log(f"  cd {Path.cwd()} && python3 master_batch_analysis.py")
    _log(f"  (Input files must be in: {in_path})")
    _log("")
    _log(f"Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    _log("=" * 70)
    print(f"\n📝 Run log: {_log_file}")

    # ─── Print folder structure ──────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  MISSION ACCOMPLISHED!")
    print(f"{'='*60}")
    print(f"  Output root: {out_root}")
    print(f"  Structure:")
    for name, d in sorted(dirs.items(), key=lambda x: x[1].name):
        n_files = sum(1 for _ in d.rglob('*'))
        print(f"    {d.name}/ ({n_files} files)")
    print()


if __name__ == "__main__":
    main()
