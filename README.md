# 🧬 Calmodulin Variant Effect Atlas (CaM‑VEA)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Build](https://github.com/QariSb/CaM-VEA/actions/workflows/build.yml/badge.svg)](https://github.com/QariSb/CaM-VEA/actions)
[![Python](https://img.shields.io/badge/python-3.11%2B-brightgreen)]()
[![Data](https://img.shields.io/badge/data-ClinVar%20%7C%20ICalmR%20%7C%20UniProt-orange)]()

---

### 🧩 Overview
**CaM‑VEA (Calmodulin Variant Effect Atlas)** is a full‑coverage *in silico* saturation mutagenesis and clinical annotation project for the **human calmodulin protein (CALM1/2/3)**.  
It provides physicochemical features, EF‑hand context, and integrated ClinVar variant metadata for all **2,831 possible missense mutations** (149 amino acids × 19 substitutions).

---

### 🔍 Objectives
- Map **every single amino‑acid substitution** in calmodulin to EF‑hand motifs and biophysical features.  
- Merge **ClinVar + International Calmodulinopathy Registry (ICalmR)** data for variant labels.  
- Enable training and benchmarking of variant effect predictors (ΔΔG, AlphaMissense, ESM, etc.).  
- Provide an open, standardized dataset for **variant interpretation and mechanism inference**.

---

### 📂 Project structure
```
CaM-VEA/
├── README.md
├── LICENSE
├── .gitignore
├── .github/workflows/build.yml
├── data/
│   ├── calm_human.fasta
│   ├── ef_hand_annotation.csv
│   ├── amino_acid_properties.csv
│   ├── clinvar_miner_CALM1.csv
│   ├── clinvar_miner_CALM2.csv
│   └── clinvar_miner_CALM3.csv
├── scripts/
│   └── 01_build_variant_space.py
├── outputs/
│   ├── cam_saturation_variants.csv
│   ├── cam_saturation_variants_with_clinvar.csv
│   └── cam_saturation_variants_labeled.csv
└── notebooks/
```

---

### ⚙️ Reproducibility
The repository builds automatically on GitHub Actions (Python 3.11) using:

```bash
python scripts/01_build_variant_space.py
```

To rebuild locally:

```bash
pip install pandas numpy biopython
python scripts/01_build_variant_space.py
```

---

### 🧠 Data sources
- **ClinVar** (2025) – curated CALM1/2/3 pathogenic/likely‑pathogenic variants  
- **ICalmR (International Calmodulinopathy Registry)** – case series & phenotype mapping  
- **UniProt P0DP23 / P0DP24 / P0DP25** – canonical sequence references  
- **Halling et al., PNAS 2016** – EF‑hand motif positions and Ca²⁺ coordination residues  

---

### 📊 Output columns
| Column | Description |
|---------|-------------|
| `position` | Residue number (1–149) |
| `wt_aa`, `alt_aa` | Wild‑type / mutant amino acids |
| `region`, `ef_loop`, `is_ca_coord` | EF‑hand / lobe context |
| `delta_charge`, `delta_hydropathy`, `delta_volume_A3` | Physicochemical deltas |
| `mechanism_prior` | Mechanistic prior class |
| `ClinVar_Variation_ID`, `Condition`, `Review_status`, `PMIDs` | Integrated ClinVar metadata |

---

### 🧾 Citation
If you use this dataset or code, please cite:

> **Basit A. (2025).** *Calmodulin Variant Effect Atlas (CaM‑VEA): a complete saturation and ClinVar‑integrated resource for CALM1–3 variant interpretation.*  
> [https://github.com/QariSb/CaM-VEA](https://github.com/QariSb/CaM-VEA)

---

### 🙌 Acknowledgments
Project maintained by **Abdul Basit (QariSb)**.  
Developed with the assistance of GPT‑5 for reproducible computational biology workflows.
