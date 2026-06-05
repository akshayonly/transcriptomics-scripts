#!/usr/bin/env python
# coding: utf-8

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Patch

os.makedirs("results", exist_ok=True)

# =============================================================================
# PATHS — edit here if your project is in a different location
# =============================================================================
BASE        = "/teamspace/studios/this_studio/2023_paper"
COUNTS_FILE = f"{BASE}/counts/counts_raw.txt"
SRA_FILE    = f"{BASE}/SraRunTable.csv"
SAMPLE_FILE = f"{BASE}/SampleInfo.csv"

# =============================================================================
# 1. LOAD COUNT MATRIX
# =============================================================================
df     = pd.read_csv(COUNTS_FILE, sep="\t", comment="#", index_col=0)
annot  = df.iloc[:, :5].copy()   # Chr, Start, End, Strand, Length
counts = df.iloc[:, 5:].copy()

# Clean column names: extract SRR accession from full BAM path
counts.columns = [re.search(r'(SRR\d+)', c).group(1) for c in counts.columns]

print(f"Genes   : {counts.shape[0]:,}")
print(f"Samples : {counts.shape[1]}")
counts.head()

# =============================================================================
# 2. METADATA
# =============================================================================
sra_info    = pd.read_csv(SRA_FILE)
sample_info = pd.read_csv(SAMPLE_FILE)

sample_sra  = sra_info[["Run", "LibraryName"]]
sample_info.rename(columns={"gsm_accession": "LibraryName",
                             "sample_title":  "Sample"}, inplace=True)
metadata = pd.merge(sample_info, sample_sra)

def parse_condition(title):
    if title.count("_") == 2:          # BOP27_ndh_1 or BOP27_ndh_2
        return "parent"
    for lin in ["A13", "A14", "A15", "A16"]:
        if lin in title:
            return f"evolved_{lin}"
    return "unknown"

def parse_lineage(title):
    for lin in ["A13", "A14", "A15", "A16"]:
        if lin in title:
            return lin
    return "parent"

def parse_rep(title):
    return int(title[-1])

def make_label(row):
    return f"{row['lineage']}_rep{row['rep']}"

metadata["condition"] = metadata["Sample"].apply(parse_condition)
metadata["lineage"]   = metadata["Sample"].apply(parse_lineage)
metadata["rep"]       = metadata["Sample"].apply(parse_rep)
metadata["label"]     = metadata.apply(make_label, axis=1)
metadata = metadata.set_index("Run")

# Reorder counts to match metadata
counts       = counts[metadata.index]
srr_to_label = metadata["label"].to_dict()

metadata[["LibraryName", "Sample", "condition", "lineage", "rep", "label"]]

# =============================================================================
# 3. COLOUR PALETTE
# =============================================================================
PALETTE = {
    "parent":       "#984ea3",
    "evolved_A13":  "#e41a1c",
    "evolved_A14":  "#ff7f00",
    "evolved_A15":  "#4daf4a",
    "evolved_A16":  "#377eb8",
}

# =============================================================================
# 4. CPM NORMALISATION
# =============================================================================
lib_sizes = counts.sum()

# CPM — corrects for library size
cpm     = counts.div(lib_sizes / 1e6)

# log2 CPM with pseudocount
log_cpm = np.log2(cpm + 1)

# Filter: keep genes with CPM > 1 in at least 2 samples
keep         = (cpm > 1).sum(axis=1) >= 2
log_cpm_filt = log_cpm[keep]

# Display-labelled versions (SRR → A13_rep1 etc.)
counts_display   = counts.rename(columns=srr_to_label)
log_cpm_display  = log_cpm_filt.rename(columns=srr_to_label)

print(f"Genes after CPM>1 filter : {log_cpm_filt.shape[0]:,} / {log_cpm.shape[0]:,}")
print(f"\nLibrary sizes:")
for srr, n in lib_sizes.items():
    print(f"  {metadata.loc[srr,'label']:<20} {n:>12,.0f}   ({metadata.loc[srr,'condition']})")

# Save normalised counts with readable column names
cpm_out = cpm.rename(columns=srr_to_label).round(3)
cpm_out.to_csv("results/py_counts_normalized.csv")
print("\nSaved: results/py_counts_normalized.csv")

# =============================================================================
# 5. LIBRARY SIZE BAR CHART
# =============================================================================
lib_sizes_display = counts_display.sum()
colors = [PALETTE[metadata.loc[srr, "condition"]] for srr in counts.columns]

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(lib_sizes_display.index, lib_sizes_display.values / 1e6,
       color=colors, edgecolor="white", width=0.7)

ax.set_ylabel("Assigned read pairs (millions)")
ax.set_title("Library sizes per sample", fontweight="bold")
ax.set_xticklabels(lib_sizes_display.index, rotation=45, ha="right", fontsize=9)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.1f} M"))

legend_elements = [Patch(facecolor=v, label=k) for k, v in PALETTE.items()]
ax.legend(handles=legend_elements, title="Condition",
          bbox_to_anchor=(1, 1), loc="upper left", fontsize=9)

plt.tight_layout()
plt.savefig("results/py_libsize.pdf", bbox_inches="tight")
plt.show()

# =============================================================================
# 6. PCA
# =============================================================================
scaler = StandardScaler()
X      = scaler.fit_transform(log_cpm_display.T)   # samples × genes

pca     = PCA(n_components=min(5, X.shape[0]))
coords  = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(
    coords[:, :2],
    columns=["PC1", "PC2"],
    index=log_cpm_display.columns          # readable labels
)
pca_df["condition"] = metadata["condition"].values
pca_df["lineage"]   = metadata["lineage"].values
pca_df["rep"]       = metadata["rep"].values

fig, ax = plt.subplots(figsize=(7, 5))
for cond, grp in pca_df.groupby("condition"):
    ax.scatter(grp["PC1"], grp["PC2"],
               color=PALETTE[cond], s=140, zorder=3, label=cond,
               edgecolors="white", linewidths=0.8)
    for label, row in grp.iterrows():
        ax.annotate(label,                     # e.g. "A13_rep1"
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
plt.show()

# =============================================================================
# 7. SAMPLE-SAMPLE SPEARMAN CORRELATION HEATMAP
# =============================================================================
corr_matrix = log_cpm_display.corr(method="spearman")

fig, ax = plt.subplots(figsize=(9, 8))
sns.heatmap(
    corr_matrix,
    annot=True, fmt=".3f", annot_kws={"size": 8},
    cmap="RdYlBu_r", vmin=0.85, vmax=1.0,
    linewidths=0.5, linecolor="white",
    ax=ax
)
ax.set_title("Sample–sample Spearman correlation\n(log2 CPM, filtered genes)",
             fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)

plt.tight_layout()
plt.savefig("results/py_correlation.pdf", bbox_inches="tight")
plt.show()

# =============================================================================
# 8. TOP 50 MOST VARIABLE GENES HEATMAP
# =============================================================================
gene_var = log_cpm_display.var(axis=1)
top50    = gene_var.nlargest(50).index
mat      = log_cpm_display.loc[top50]

# Row z-score
mat_z = mat.subtract(mat.mean(axis=1), axis=0).divide(mat.std(axis=1), axis=0)

# Column colour annotation strip
col_colors_list = [PALETTE[metadata.loc[srr, "condition"]] for srr in counts.columns]
col_colors_ser  = pd.Series(col_colors_list, index=log_cpm_display.columns)

fig, ax = plt.subplots(figsize=(10, 14))
sns.heatmap(
    mat_z,
    cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
    linewidths=0.3, linecolor="grey",
    xticklabels=True, yticklabels=True,
    ax=ax
)
ax.set_title("Top 50 most variable genes\n(log2 CPM, row z-scored)",
             fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

plt.tight_layout()
plt.savefig("results/py_heatmap_top50.pdf", bbox_inches="tight")
plt.show()

# =============================================================================
# 9. SANITY CHECK — ndh counts (deleted gene, should be ~0)
# =============================================================================
if "ndh" in counts.index:
    ndh_counts = counts.rename(columns=srr_to_label).loc["ndh"]
    print("ndh raw counts (gene is deleted — expect ~0):")
    print(ndh_counts.to_string())
else:
    print("'ndh' not found in count matrix")
