#!/usr/bin/env python3
"""
Full Self-contained Script
--------------------------
1. Reads CALM_all_genes_residue_with_benign.csv and aligned_sequences_trimmed.fasta
2. Computes Shannon entropy per alignment column
3. Maps entropy to residue positions
4. Computes normalized conservation = 1 - H/log2(20)
5. Applies theoretical cutoff (0.5)
6. Plots conservation vs residue position and labels misclassified residues
7. Exports the misclassification table and the figure
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import math
import os

# ================================
# 1️⃣ File paths
# ================================
CSV_PATH = "CALM_all_genes_residue_with_benign.csv"
FASTA_PATH = "data/aligned_sequences_trimmed.fasta"
OUT_PNG = "conservation_vs_position_misclassified.png"
OUT_CSV = "misclassified_residues_by_conservation_cutoff.csv"

# ================================
# 2️⃣ FASTA parser
# ================================
def read_fasta(path):
    records = []
    with open(path, "r") as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records

records = read_fasta(FASTA_PATH)
fa_df = pd.DataFrame(records, columns=["id", "sequence"])
fa_df["length"] = fa_df["sequence"].str.len()

# Verify alignment is uniform
unique_lengths = sorted(fa_df["length"].unique())
if len(unique_lengths) != 1:
    raise ValueError(f"FASTA alignment not uniform: lengths = {unique_lengths}")
aln_len = unique_lengths[0]
n_seq = len(fa_df)

# ================================
# 3️⃣ Compute per-column entropy
# ================================
aa_list = list("ACDEFGHIKLMNPQRSTVWY")
per_pos_entropy = []
for i in range(aln_len):
    col_chars = [seq[i] for seq in fa_df["sequence"]]
    counts = Counter([c for c in col_chars if c in aa_list])
    total = sum(counts.values())
    H = 0.0
    for c in counts.values():
        p = c / total
        H -= p * math.log(p, 2)
    per_pos_entropy.append(H)

entropy_df = pd.DataFrame({
    "aln_col_1based": np.arange(1, aln_len + 1),
    "entropy_bits": per_pos_entropy
})

# ================================
# 4️⃣ Map protein residue positions to alignment columns
# ================================
ref_seq = fa_df.iloc[0]["sequence"]
prot_to_aln = {}
prot_pos = 0
for i, c in enumerate(ref_seq):
    if c != "-":
        prot_pos += 1
        prot_to_aln[prot_pos] = i + 1

# ================================
# 5️⃣ Read mutation CSV and merge entropy
# ================================
df = pd.read_csv(CSV_PATH)
required_cols = {"Gene", "Wild_residue", "Position", "Mutant_residue", "Phenotype"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in CSV: {missing}")

df = df.dropna(subset=["Position"]).copy()
df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["aln_col_1based"] = df["Position"].map(prot_to_aln)
merged = df.merge(entropy_df, on="aln_col_1based", how="left")

# Drop rows without valid entropy
merged = merged[np.isfinite(merged["entropy_bits"])].copy()

# ================================
# 6️⃣ Compute normalized conservation & apply theoretical cutoff
# ================================
merged["conservation"] = 1 - (merged["entropy_bits"] / np.log2(20))
merged["phen_lower"] = merged["Phenotype"].str.lower().str.strip()
merged = merged[merged["phen_lower"].isin(["benign", "deleterious"])].copy()

merged["pred_theoretical"] = np.where(merged["conservation"] > 0.5, "deleterious", "benign")
merged["misclassified"] = merged["pred_theoretical"] != merged["phen_lower"]
mis_df = merged[merged["misclassified"]].copy()

# ================================
# 7️⃣ Plot conservation vs position
# ================================
colors = {"benign": "tab:blue", "deleterious": "tab:red"}
plt.figure(figsize=(12, 6))
for phen, g in merged.groupby("phen_lower"):
    plt.scatter(
        g["Position"], g["conservation"], s=80, alpha=0.85,
        color=colors.get(phen, "gray"), edgecolor="black", linewidth=0.5,
        label=phen.capitalize()
    )

plt.axhline(y=0.5, color="black", lw=2, ls="--", label="Theoretical cutoff = 0.5")

# Label all misclassified residues
for _, row in mis_df.iterrows():
    if np.isfinite(row["Position"]) and np.isfinite(row["conservation"]):
        plt.text(
            row["Position"], row["conservation"] + 0.03, str(int(row["Position"])),
            fontsize=9, ha="center", va="bottom", color="red", fontweight="bold"
        )

plt.title("Normalized Conservation vs. Residue Position",
          fontsize=14, weight="bold")
plt.xlabel("Residue Position")
plt.ylim(0,1)
plt.ylabel("Normalized Conservation (1 - H/log₂(20))")
plt.legend(frameon=True, fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=600)
plt.close()

# ================================
# 8️⃣ Export misclassified residues table
# ================================
mis_df[["Gene", "Position", "Phenotype", "conservation", "pred_theoretical"]].to_csv(OUT_CSV, index=False)

print("✅ Analysis complete!")
print(f"🖼️  Plot saved as: {os.path.abspath(OUT_PNG)}")
print(f"📄 Misclassified table saved as: {os.path.abspath(OUT_CSV)}")
