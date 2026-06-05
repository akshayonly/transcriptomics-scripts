#!/usr/bin/env python3
"""
11_analyze_counts.py
Usage:   python scripts/11_analyze_counts.py
Input:   counts/counts_raw.txt
Output:  results/py_pca.pdf
         results/py_heatmap_top50.pdf
         results/py_correlation.pdf
         results/py_libsize.pdf
         results/py_counts_normalized.csv
         results/py_summary.txt

Requires: pandas, numpy, matplotlib, seaborn, scipy, sklearn
Install:  pip install pandas numpy matplotlib seaborn scipy scikit-learn
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

os.makedirs("results", exist_ok=True)

# ── Colour palette — one per condition ───────────────────────────────────────
PALETTE = {
    "parent":    "#984ea3",
    "evolved_A": "#e41a1c",
    "evolved_B": "#ff7f00",
    "evolved_C": "#4daf4a",
    "evolved_D": "#377eb8",
}

# ── 1. Load count matrix ──────────────────────────────────────────────────────
print("Loading count matrix...")
df = pd.read_csv("counts/counts_raw.txt", sep="\t", comment="#", index_col=0)

# Annotation columns: Chr, Start, End, Strand, Length
annot  = df.iloc[:, :5].copy()
counts = df.iloc[:, 5:].copy()

# Clean column names: extract SRR accession from full BAM path
# e.g. "alignments/SRR26027886.sorted.bam" → "SRR26027886"
counts.columns = [re.search(r'(SRR\d+)', c).group(1) for c in counts.columns]

print(f"  Genes   : {counts.shape[0]:,}")
print(f"  Samples : {counts.shape[1]}")
print(f"  Columns : {list(counts.columns)}\n")

# ── 2. Sample metadata ────────────────────────────────────────────────────────
metadata = pd.DataFrame({
    "sample": [
        "SRR26027894", "SRR26027895",   # parent (unevolved Δndh)
        "SRR26027892", "SRR26027893",   # evolved lineage A
        "SRR26027890", "SRR26027891",   # evolved lineage B
        "SRR26027888", "SRR26027889",   # evolved lineage C
        "SRR26027886", "SRR26027887",   # evolved lineage D
    ],
    "condition": [
        "parent",    "parent",
        "evolved_A", "evolved_A",
        "evolved_B", "evolved_B",
        "evolved_C", "evolved_C",
        "evolved_D", "evolved_D",
    ],
    "lineage": ["parent","parent","A","A","B","B","C","C","D","D"],
    "rep":     [1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
}).set_index("sample")

# Reorder counts columns to match metadata
counts = counts[metadata.index]

# ── 3. Basic summary ──────────────────────────────────────────────────────────
lib_sizes = counts.sum()
print("Library sizes (assigned read pairs):")
for s, n in lib_sizes.items():
    print(f"  {s}  {n:>12,.0f}  ({metadata.loc[s,'condition']})")

# ── 4. CPM normalisation ──────────────────────────────────────────────────────
# Counts Per Million — corrects for library size differences
cpm = counts.div(lib_sizes / 1e6)

# Log2 CPM (pseudocount +1 to handle zeros)
log_cpm = np.log2(cpm + 1)

# Filter lowly expressed genes: keep genes with CPM > 1 in at least 2 samples
keep = (cpm > 1).sum(axis=1) >= 2
log_cpm_filt = log_cpm[keep]
print(f"\nGenes after CPM>1 filter: {log_cpm_filt.shape[0]:,} / {log_cpm.shape[0]:,}\n")

# Save normalised counts
cpm.round(3).to_csv("results/py_counts_normalized.csv")
print("Saved: results/py_counts_normalized.csv")

# ── 5. Library size bar chart ─────────────────────────────────────────────────
print("Plotting library sizes...")
fig, ax = plt.subplots(figsize=(10, 4))
colors = [PALETTE[metadata.loc[s, "condition"]] for s in lib_sizes.index]
bars = ax.bar(lib_sizes.index, lib_sizes.values / 1e6, color=colors, edgecolor="white")
ax.set_ylabel("Assigned read pairs (millions)")
ax.set_title("Library sizes per sample", fontweight="bold")
ax.set_xticklabels(lib_sizes.index, rotation=45, ha="right", fontsize=9)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f M"))

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=v, label=k) for k, v in PALETTE.items()]
ax.legend(handles=legend_elements, title="Condition", bbox_to_anchor=(1, 1), loc="upper left")

plt.tight_layout()
plt.savefig("results/py_libsize.pdf", bbox_inches="tight")
plt.close()
print("Saved: results/py_libsize.pdf")

# ── 6. PCA ────────────────────────────────────────────────────────────────────
print("Plotting PCA...")
scaler = StandardScaler()
X = scaler.fit_transform(log_cpm_filt.T)   # samples × genes

pca = PCA(n_components=min(5, X.shape[0]))
coords = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(coords[:, :2], columns=["PC1", "PC2"], index=counts.columns)
pca_df["condition"] = metadata["condition"]
pca_df["lineage"]   = metadata["lineage"]
pca_df["rep"]       = metadata["rep"]

fig, ax = plt.subplots(figsize=(7, 5))
for cond, grp in pca_df.groupby("condition"):
    ax.scatter(grp["PC1"], grp["PC2"],
               color=PALETTE[cond], s=120, zorder=3, label=cond,
               edgecolors="white", linewidths=0.8)
    for _, row in grp.iterrows():
        ax.annotate(f"{row['lineage']}.{row['rep']}",
                    (row["PC1"], row["PC2"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)

ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% variance)")
ax.set_title("PCA — E. coli ndh ALE RNA-Seq\n(log2 CPM, filtered genes)",
             fontweight="bold")
ax.legend(title="Condition", bbox_to_anchor=(1, 1), loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig("results/py_pca.pdf", bbox_inches="tight")
plt.close()
print("Saved: results/py_pca.pdf")

# ── 7. Sample-sample Spearman correlation heatmap ────────────────────────────
print("Plotting correlation heatmap...")
corr_matrix = log_cpm_filt.corr(method="spearman")

# Annotation for columns/rows
col_colors = pd.Series(
    [PALETTE[metadata.loc[s, "condition"]] for s in corr_matrix.columns],
    index=corr_matrix.columns
)

fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(
    corr_matrix,
    annot=True, fmt=".3f", annot_kws={"size": 7},
    cmap="RdYlBu_r", vmin=0.9, vmax=1.0,
    linewidths=0.5, linecolor="white",
    ax=ax
)
ax.set_title("Sample–sample Spearman correlation\n(log2 CPM, filtered genes)",
             fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
plt.tight_layout()
plt.savefig("results/py_correlation.pdf", bbox_inches="tight")
plt.close()
print("Saved: results/py_correlation.pdf")

# ── 8. Heatmap — top 50 most variable genes ──────────────────────────────────
print("Plotting top-50 variable genes heatmap...")
gene_var   = log_cpm_filt.var(axis=1)
top50      = gene_var.nlargest(50).index
mat        = log_cpm_filt.loc[top50]

# Z-score across samples (row-scale)
mat_z = mat.subtract(mat.mean(axis=1), axis=0).divide(mat.std(axis=1), axis=0)

# Column colour bar
col_colors_list = [PALETTE[metadata.loc[s, "condition"]] for s in mat_z.columns]

fig, ax = plt.subplots(figsize=(10, 14))
sns.heatmap(
    mat_z,
    cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
    linewidths=0.3, linecolor="grey",
    xticklabels=True, yticklabels=True,
    ax=ax
)
ax.set_title("Top 50 most variable genes\n(log2 CPM, row z-scored)", fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

# Add condition colour strip above heatmap
for i, (s, color) in enumerate(zip(mat_z.columns, col_colors_list)):
    ax.add_patch(plt.Rectangle((i, len(top50)), 1, 1.5,
                                color=color, clip_on=False,
                                transform=ax.get_transform()))

plt.tight_layout()
plt.savefig("results/py_heatmap_top50.pdf", bbox_inches="tight")
plt.close()
print("Saved: results/py_heatmap_top50.pdf")

# ── 9. Sanity check — ndh counts ─────────────────────────────────────────────
print("\nSanity check — ndh gene counts (should be ~0, gene is deleted):")
if "ndh" in counts.index:
    print(counts.loc["ndh"].to_string())
else:
    print("  'ndh' not found in count matrix index")

# ── 10. Summary report ────────────────────────────────────────────────────────
summary_lines = [
    "=== counts_raw.txt Python Analysis Summary ===\n",
    f"Total genes in matrix  : {counts.shape[0]:,}",
    f"Samples                : {counts.shape[1]}",
    f"Genes after CPM filter : {log_cpm_filt.shape[0]:,}",
    "",
    "Library sizes:",
]
for s, n in lib_sizes.items():
    summary_lines.append(f"  {s}  {n:>12,.0f}  {metadata.loc[s,'condition']}")

summary_lines += [
    "",
    "PCA variance explained:",
    f"  PC1: {var_exp[0]:.1f}%",
    f"  PC2: {var_exp[1]:.1f}%",
    f"  PC3: {var_exp[2]:.1f}%",
    "",
    "Outputs saved to results/:",
    "  py_libsize.pdf           — library size bar chart",
    "  py_pca.pdf               — PCA plot",
    "  py_correlation.pdf       — sample-sample Spearman correlation",
    "  py_heatmap_top50.pdf     — top 50 variable genes heatmap",
    "  py_counts_normalized.csv — CPM normalised counts",
]

summary = "\n".join(summary_lines)
print("\n" + summary)
with open("results/py_summary.txt", "w") as f:
    f.write(summary)
print("\nSaved: results/py_summary.txt")
print("\nDone.")
