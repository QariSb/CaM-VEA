#!/usr/bin/env python3
"""
Background-Weighted Conservation with Henikoff Sequence Weighting
-----------------------------------------------------------------
1. Reads CALM_all_genes_residue_with_benign.csv and aligned_sequences_trimmed.fasta
2. Applies Henikoff sequence weighting to downweight redundant sequences
3. Computes background-weighted relative entropy (IC)
4. Normalizes conservation: IC_norm = IC / log2(1/q_min)
5. Maps entropy to residue positions
6. Applies theoretical cutoff (0.5)
7. Plots conservation vs residue position and labels misclassified residues
8. Exports misclassification table and figure

Author: Abdul Basit + ChatGPT (2025)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import math, os

# ================================
# 1️⃣ File paths
# ================================
CSV_PATH = "CALM_all_genes_residue_with_benign.csv"
FASTA_PATH = "data/aligned_sequences_trimmed.fasta"
OUT_PNG = "conservation_bg_henikoff_vs_position_misclassified.png"
OUT_CSV = "misclassified_residues_bg_henikoff.csv"

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

unique_lengths = sorted(fa_df["length"].unique())
if len(unique_lengths) != 1:
    raise ValueError(f"FASTA alignment not uniform: lengths = {unique_lengths}")
aln_len = unique_lengths[0]
n_seq = len(fa_df)

print(f"✅ Loaded alignment with {n_seq} sequences, length = {aln_len}")

# ================================
# 3️⃣ Henikoff sequence weighting
# ================================
aa_list = list("ACDEFGHIKLMNPQRSTVWY")
AA2IDX = {a:i for i,a in enumerate(aa_list)}

def henikoff_weights(seqs):
    """
    Compute Henikoff sequence weights for an aligned set of sequences.
    Returns array of weights (len = n_seq), normalized to sum to 1.
    """
    n_seq = len(seqs)
    L = len(seqs[0])
    weights = np.zeros(n_seq, dtype=float)

    # Iterate over alignment columns
    for j in range(L):
        col = [s[j] for s in seqs]
        residues = [c for c in col if c != "-"]
        if not residues:
            continue
        unique_res = sorted(set(residues))
        r = len(unique_res)
        counts = Counter(residues)
        for i, s in enumerate(seqs):
            c = s[j]
            if c == "-" or c not in counts:
                continue
            weights[i] += 1.0 / (r * counts[c])

    # Normalize so that total weight = number of sequences
    weights = weights / np.sum(weights) * n_seq
    weights /= np.sum(weights)  # ensure normalized to 1
    return weights

seqs = fa_df["sequence"].tolist()
seq_weights = henikoff_weights(seqs)
print(f"✅ Henikoff weights computed (sum={seq_weights.sum():.3f})")

# ================================
# 4️⃣ Background-weighted conservation (relative entropy)
# ================================
def compute_background_q_weighted(seqs, weights, alpha=1.0):
    """Weighted background frequencies with Laplace smoothing."""
    counts = np.zeros(len(aa_list), dtype=float)
    for w, s in zip(weights, seqs):
        for c in s:
            if c in AA2IDX:
                counts[AA2IDX[c]] += w
    total = counts.sum()
    q = (counts + alpha) / (total + alpha * len(aa_list))
    return q

def column_relative_entropy_weighted(seqs, weights, q, alpha=0.5):
    """
    Weighted per-column IC = sum_a p(a) * log2(p(a)/q(a)).
    Returns: (ic_bits, ic_norm)
    """
    L = len(seqs[0])
    ic_bits = np.zeros(L, dtype=float)
    q = np.asarray(q, dtype=float)
    q_min = float(np.min(q))
    ic_max = math.log2(1.0 / q_min)

    for i in range(L):
        col_counts = np.zeros(len(aa_list), dtype=float)
        for w, s in zip(weights, seqs):
            c = s[i]
            if c in AA2IDX:
                col_counts[AA2IDX[c]] += w
        total = col_counts.sum()
        if total == 0:
            ic_bits[i] = 0.0
            continue
        p = (col_counts + alpha) / (total + alpha * len(aa_list))
        ic = float(np.sum(p * (np.log2(p) - np.log2(q))))
        ic_bits[i] = max(ic, 0.0)

    ic_norm = ic_bits / ic_max if ic_max > 0 else np.zeros_like(ic_bits)
    return ic_bits, ic_norm

# Weighted background
q_empirical = compute_background_q_weighted(seqs, seq_weights, alpha=1.0)
ic_bits, ic_norm = column_relative_entropy_weighted(seqs, seq_weights, q_empirical, alpha=0.5)

entropy_df = pd.DataFrame({
    "aln_col_1based": np.arange(1, len(ic_bits) + 1),
    "ic_bits": ic_bits,
    "conservation_bg": ic_norm
})

# ================================
# 5️⃣ Map protein residue positions to alignment columns
# ================================
ref_seq = fa_df.iloc[0]["sequence"]
prot_to_aln = {}
prot_pos = 0
for i, c in enumerate(ref_seq):
    if c != "-":
        prot_pos += 1
        prot_to_aln[prot_pos] = i + 1

# ================================
# 6️⃣ Read mutation CSV and merge entropy
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
merged = merged[np.isfinite(merged["conservation_bg"])].copy()

# ================================
# 7️⃣ Apply cutoff and classify
# ================================
merged["phen_lower"] = merged["Phenotype"].str.lower().str.strip()
merged = merged[merged["phen_lower"].isin(["benign", "deleterious"])].copy()
merged["pred_theoretical"] = np.where(merged["conservation_bg"] > 0.5, "deleterious", "benign")
merged["misclassified"] = merged["pred_theoretical"] != merged["phen_lower"]
mis_df = merged[merged["misclassified"]].copy()

n_total = len(merged)
n_mis = len(mis_df)
print(f"🔹 Total residues classified: {n_total}")
print(f"🔹 Misclassified: {n_mis} ({100*n_mis/n_total:.1f}%)")

# ================================
# 8️⃣ Plot conservation vs position
# ================================
colors = {"benign": "tab:blue", "deleterious": "tab:red"}
plt.figure(figsize=(12, 6))
for phen, g in merged.groupby("phen_lower"):
    plt.scatter(
        g["Position"], g["conservation_bg"], s=80, alpha=0.85,
        color=colors.get(phen, "gray"), edgecolors="black", linewidth=0.5,
        label=phen.capitalize()
    )

plt.axhline(y=0.5, color="black", lw=2, ls="--", label="Theoretical cutoff = 0.5")
for _, row in mis_df.iterrows():
    if np.isfinite(row["Position"]) and np.isfinite(row["conservation_bg"]):
        plt.text(
            row["Position"], row["conservation_bg"] + 0.03, str(int(row["Position"])),
            fontsize=9, ha="center", va="bottom", color="red", fontweight="bold"
        )

plt.title("Henikoff-Weighted Background Conservation vs. Residue Position",
          fontsize=14, weight="bold")
plt.xlabel("Residue Position")
plt.ylabel("Normalized Conservation (IC / log₂(1/q_min))")
plt.ylim(0,1)
plt.legend(frameon=True, fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=600)
plt.close()

# ================================
# 9️⃣ Export misclassified residues table
# ================================
mis_df_out = mis_df.copy()
mis_df_out["conservation_bg"] = mis_df_out["conservation_bg"].round(3)
mis_df_out["ic_bits"] = mis_df_out["ic_bits"].round(3)
mis_df_out[["Gene", "Position", "Phenotype", "conservation_bg", "pred_theoretical"]].to_csv(OUT_CSV, index=False)

print("✅ Analysis complete!")
print(f"🖼️  Plot saved as: {os.path.abspath(OUT_PNG)}")
print(f"📄 Misclassified table saved as: {os.path.abspath(OUT_CSV)}")
