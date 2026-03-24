# Protein Property analysis

A batch analysis pipeline that evaluates protein/antibody candidates across multiple biophysical properties from PDB/CIF structure files.

## Overview

Given a set of protein structure files, this pipeline automatically runs 6 analysis tools in sequence and generates a unified report with PASS/FAIL interpretation for each candidate.

**Evaluated properties:**

* **Thermal stability** — ThermoMPNN (ΔΔG, kcal/mol)
* **Structural solubility** — GATSol (probability)
* **Aggregation propensity** — Aggrescan (Na4vSS, nHS, AAT) + AggNet (APR count/fraction, peak score)
* **Expression solubility** — NetSolP (E. coli solubility \& usability)
* **Immunogenicity** — MHC-I (HLA-A/B/C alleles) and MHC-II (HLA-DRB1 alleles) via sliding-window %Rank

## Example Output

Results are saved as Excel and TSV with inline PASS/BORDERLINE/CAUTION/FAIL interpretation per metric.

|ID|Stability\_ddG|Solubility|Aggrescan\_Na4vSS|NetSolP\_Sol|
|-|-|-|-|-|
|9I6Q|-0.007 |0.522 |-9.27 |0.358 |
|WT|-0.005 |0.487 |29.76 |0.289 |

A heatmap visualization is also generated for quick visual comparison across candidates.

## Pipeline Architecture

```
Input: .pdb / .cif files
        │
        ▼
┌─────────────────────────────┐
│  Structure Preprocessing    │  BioPython — chain parsing, renumbering, HETATM removal
└────────────┬────────────────┘
             │
     ┌───────┴────────┐
     ▼                ▼
\[1] ThermoMPNN    \[2] GATSol        Stability \& solubility (structure-based)
\[3] MHC Pipeline  \[4] Aggrescan     Immunogenicity \& aggregation (sequence-based)
\[5] AggNet        \[6] NetSolP       Aggregation (DL) \& E.coli expression
     │
     ▼
┌─────────────────────────────┐
│  Summary Report             │  Excel + TSV with interpretation + heatmaps
└─────────────────────────────┘
```

## Requirements

The following tools must be installed and configured on the system:

|Tool|Purpose|Environment|
|-|-|-|
|ThermoMPNN|Stability ΔΔG|`/root/ThermoMPNN`|
|GATSol|Structural solubility|`/root/GATSol/Predict`|
|Aggrescan|Classical aggregation|`/root/aggrescan`|
|AggNet|DL-based aggregation|`/root/AggNet`|
|NetSolP|E. coli solubility|conda env `netsolp`|
|MHC predictor|Immunogenicity|`/root/mhc1\_env`|

Python dependencies: `biopython`, `pandas`, `matplotlib`, `seaborn`, `openpyxl`

## Usage

1. Place input `.pdb` or `.cif` files in the input directory (default: `/root/For\_phyproper`)
2. Optionally name the reference/wildtype file with a tag: `WT`, `ref`, `test`, or `01`
3. Run:

```bash
python3 master\_batch\_analysis.py
```

Output is saved to a timestamped folder with the following structure:

```
output\_YYYYMMDD\_HHMMSS/
├── 01\_ThermoMPNN/
├── 02\_GATSol/
├── 03\_MHC/
├── 04\_AGGRESCAN/
├── 05\_AggNet/
├── 06\_NetSolP/
└── Summary/
    ├── Final\_Assembly\_Report.xlsx   ← Main report with interpretation
    ├── Final\_Assembly\_Report.tsv
    ├── Final\_Assembly\_Report\_raw.tsv
    └── Visualizations/              ← Heatmaps
```

## Acceptance Criteria

|Metric|Direction|Threshold|
|-|-|-|
|Stability\_ddG|lower|< 0.0 kcal/mol|
|Solubility (GATSol)|higher|> 0.5|
|Aggrescan\_Na4vSS|lower|< 0.0|
|Aggrescan\_nHS|lower|< 5|
|AggNet\_APR\_frac|lower|< 0.05|
|NetSolP\_Sol|higher|> 0.5|
|MHC-I/II avg %Rank|higher|> 40|



