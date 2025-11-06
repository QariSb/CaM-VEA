"""
01_build_variant_space.py
-------------------------

Calmodulin Variant Effect Atlas (CaM-VEA)
-----------------------------------------
Generates the complete single-missense saturation dataset for the
human calmodulin canonical sequence (CALM1/2/3; UniProt P0DP23/24/25).

Steps:
1. Load canonical calmodulin sequence (149 aa)
2. Annotate residues by EF-hand and lobe
3. Compute all 19 alternative amino-acid substitutions
4. Calculate physicochemical deltas (Δcharge, Δhydropathy, Δvolume)
5. Assign mechanistic priors (e.g. “Ca2+-coordination_site”)
6. Save to outputs/cam_saturation_variants.csv
"""

import pandas as pd
from pathlib import Path

# -------------------- Paths --------------------
root = Path(__file__).resolve().parents[1]
data_dir = root / "data"
out_dir = root / "outputs"
out_dir.mkdir(parents=True, exist_ok=True)

fasta_path = data_dir / "calm_human.fasta"
ef_path = data_dir / "ef_hand_annotation.csv"
aa_path = data_dir / "amino_acid_properties.csv"
out_path = out_dir / "cam_saturation_variants.csv"

# -------------------- Load inputs --------------------
ef_df = pd.read_csv(ef_path)
aa_df = pd.read_csv(aa_path)
aa_map = aa_df.set_index("aa").to_dict(orient="index")

# Load WT sequence
lines = [l.strip() for l in open(fasta_path) if l.strip() and not l.startswith(">")]
seq = "".join(lines)
assert len(seq) == 149, f"Expected 149 aa, got {len(seq)}"
assert all(c in aa_map for c in seq), "Invalid amino acid in sequence"

# -------------------- Build variant space --------------------
aas = list(aa_df["aa"].values)
variants = []
for pos in range(1, 150):
    wt = seq[pos - 1]
    for alt in aas:
        if alt == wt:
            continue
        variants.append({"position": pos, "wt_aa": wt, "alt_aa": alt})

df = pd.DataFrame(variants).merge(ef_df, on="position", how="left")

# -------------------- Compute physicochemical deltas --------------------
def delta(a, b, field):
    return aa_map[b][field] - aa_map[a][field]

for field in ["charge", "hydropathy", "volume_A3"]:
    df[f"delta_{field}"] = df.apply(lambda r: delta(r["wt_aa"], r["alt_aa"], field), axis=1)

# -------------------- Mechanistic context heuristic --------------------
def mech_tag(r):
    if r["ef_loop"] and r["is_ca_coord"]:
        return "Ca2+-coordination_site"
    if r["ef_loop"] and r["is_conserved_gly"]:
        return "EF-loop_flex_gly"
    if r["ef_loop"]:
        return "EF-loop_other"
    if r["region"] == "Linker":
        return "Interlobe_linker"
    return f"{r['region']}_surface"

df["mechanism_prior"] = df.apply(mech_tag, axis=1)

# -------------------- Save output --------------------
df.to_csv(out_path, index=False)
print(f"[✅] Saved {len(df):,} variants to {out_path}")
print(df.head(10))
