#imports
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pysam
from itertools import chain
import argparse

#parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument("--summary_folder", required=True)
parser.add_argument("--snv_file", required=True)
parser.add_argument("--genes_coordinates_file", required=True)
parser.add_argument("--sample_name", required=True)
args = parser.parse_args()

#getting input tables
genes_table = pd.read_csv(args.genes_coordinates_file, sep = '\t', header = None)
genes_table.columns = ['chr','start','end','gene_symbol','ensid']
genes_table["chr_num"] = genes_table["chr"].str.replace("chr", "").astype(int)
genes_table = genes_table.sort_values(by=['chr_num','start'])
genes_table = genes_table.loc[:,['chr','start','end','gene_symbol','ensid']].reset_index()
genes_table = genes_table[genes_table['gene_symbol'] != 'EGFR-AS1']

CNV_summary = pd.read_csv(os.path.join(args.summary_folder,'CNV_overlapping_genes.tsv'), sep = '\t', header = None)
CNV_summary = CNV_summary.iloc[:,[7,8,14,16]]
CNV_summary.columns = ['Tot_CN', 'Minor_CN', 'Gene', 'Overlap(bp)']
CNV_summary = CNV_summary[CNV_summary['Gene'] != 'EGFR-AS1']

SV_summary = pd.read_csv(os.path.join(args.summary_folder,'SV_overlapping_genes.vcf'), sep = '\t', header = None)

genes_SV_overlapped = SV_summary.iloc[:,13]
unique_genes = list(set(genes_SV_overlapped))
counts_SV_by_gene = [genes_SV_overlapped.tolist().count(gene) for gene in unique_genes]
SV_counts = dict(zip(unique_genes,counts_SV_by_gene))
SV_counts = [SV_counts[gene] if gene in SV_counts.keys() else 0 for gene in genes_table['gene_symbol']]
SV_counts = pd.DataFrame({'genes':genes_table['gene_symbol'], 'counts': SV_counts})
SV_counts = SV_counts[SV_counts['genes'] != 'EGFR-AS1']
SV_counts = SV_counts.sort_values(by = 'genes')

snp_summary_path = args.snv_file
vcf = pysam.VariantFile(snp_summary_path)


##formatting input tables##
#snv#
desc = vcf.header.info["CSQ"].description
fields = desc.split("Format: ")[-1].split("|")
gene_idx = fields.index("SYMBOL")   
effect_idx = fields.index("Consequence")

genes = list()
effects = list()

for rec in vcf:
    for csq in rec.info.get("CSQ", []):
        parts = csq.split("|")
        if len(parts) > gene_idx and parts[gene_idx]:
            genes.append(parts[gene_idx])
            effects.append(parts[effect_idx])

hit = [i for i, g in enumerate(genes) if g in genes_table['gene_symbol'].tolist()]
snp_df = pd.DataFrame({'gene': [genes[h] for h in hit], 'effect': [effects[h] for h in hit]})

possible_effects = list(set(snp_df['effect']))

effects_vec = list()
effects_counts = list()
genes_vec = list()

for g in genes_table['gene_symbol']:
    genes_vec.append([g]*len(possible_effects))
    effects_vec.append(possible_effects)
    eff_current_counts = [0]*len(possible_effects)
    if g in snp_df['gene'].tolist() :
        snp_df_tmp = snp_df[snp_df['gene']==g]
        for eff in snp_df_tmp['effect']:
            eff_current_counts[possible_effects.index(eff)] += 1
    effects_counts.append(eff_current_counts)

genes_vec = list(chain.from_iterable(genes_vec))
effects_vec = list(chain.from_iterable(effects_vec))
effects_counts = list(chain.from_iterable(effects_counts))

snp_df_forplot = pd.DataFrame({'gene': genes_vec, 'effects': effects_vec, 'counts':effects_counts})
snp_df_forplot = snp_df_forplot.sort_values(by = 'gene')
#CNV#
genes_length = (
    genes_table
    .assign(length=genes_table["end"] - genes_table["start"])
    .set_index("gene_symbol")["length"]
)

CNV_summary["Overlap(bp)"] = (
    CNV_summary["Overlap(bp)"] /
    CNV_summary["Gene"].map(genes_length)
)

CNV_summary["Tot_CN"] = CNV_summary["Tot_CN"].round().astype("Int64")
CNV_summary["Minor_CN"] = CNV_summary["Minor_CN"].round().astype("Int64")
CNV_summary.columns = ['Tot_CN', 'Minor_CN', 'Gene','Normalized overlap']

def handle_duplicated_genes_spaces(df, col="Gene"):
    df = df.copy()
    seen = {}  

    out = []
    for g in df[col].astype(str):
        n = seen.get(g, 0)
        out.append(g + (" " * n))
        seen[g] = n + 1

    df[col] = out
    return df

CNV_summary = handle_duplicated_genes_spaces(CNV_summary)
CNV_summary_to_plot = CNV_summary.melt(id_vars='Gene', var_name='feature', value_name='features_values')
CNV_summary_to_plot = CNV_summary_to_plot.sort_values(by = 'Gene')
##plotting##
sns.set_theme(
    context="paper",
    style="whitegrid",
    font="sans-serif",
    font_scale=1.1
)
plt.rcParams["figure.dpi"] = 100
plt.rcParams.update({
    "grid.linestyle": "--",
    "grid.alpha": 0.3,
    "axes.linewidth": 1.0,
    "axes.edgecolor": "black",
})
#SV#
plt.figure(figsize=(15, 5))
plt.tick_params(
    axis="x",
    which="major",
    length=8,
    width=1.2,
    bottom=True,
    top=False
)
sns.barplot(data=SV_counts, x="genes", y="counts",palette="muted",hue = 'counts',legend = False, edgecolor = 'black')
plt.xticks(rotation=60, fontsize = 9, fontweight = 'bold')
plt.margins(x=0.025)
plt.tight_layout()
plt.xlabel("Gene", fontsize=10)
plt.ylabel("SV count", fontsize=10)
plt.title(args.sample_name + " SVs count per gene", fontsize=14, fontstyle = 'italic')

plt.savefig(args.summary_folder + "SV_plot.pdf", bbox_inches="tight")

#snv#
plt.figure(figsize=(20, 5), dpi = 100)
sns.barplot(data=snp_df_forplot, x="gene", y="counts", hue="effects", width=1, edgecolor = 'black')
plt.legend(fontsize=8)
plt.tight_layout()
plt.xlabel("Gene", fontsize=10)
plt.ylabel("SNV count", fontsize=10)
plt.title(args.sample_name + " snv count and effect", fontsize=14, fontstyle = 'italic')
plt.xticks(rotation=60, fontsize = 9, fontweight = 'bold')
plt.margins(x=0.02)
plt.tick_params(
    axis="x",
    which="major",
    length=8,
    width=1.2,
    bottom=True,
    top=False
)

plt.savefig(args.summary_folder + "snv_plot.pdf", bbox_inches="tight")
#cnv#
FIGSIZE = (30, 8)
DPI = 100
XROT = 60
XTICK_FONTSIZE = 9
XLABEL_FONTSIZE = 12
TITLE_FONTSIZE = 17
XMARGIN = 0.03

fig, ax1 = plt.subplots(figsize=FIGSIZE, dpi=DPI)
CNV_summary_to_plot = CNV_summary_to_plot.sort_values(by = 'Gene')

plt.tick_params(
    axis="x",
    which="major",
    length=8,
    width=1.2,
    bottom=True,
    top=False
)

sns.barplot(
    data=CNV_summary_to_plot[CNV_summary_to_plot["feature"].isin(["Tot_CN", "Minor_CN"])],
    x="Gene",
    y="features_values",
    hue="feature",
    width=0.4,
    edgecolor="black",
    errorbar=None,
    ax=ax1
)

ax1.set_ylabel("Copy number", color = 'darkorange')
ax1.tick_params(axis='y', colors = 'blue')


ax2 = ax1.twinx()
color = "#55A868"

sns.barplot(
    data=CNV_summary_to_plot[CNV_summary_to_plot["feature"] == "Normalized overlap"],
    x="Gene",
    y="features_values",
    color=color,
    width=0.2,
    edgecolor="black",
    errorbar=None,
    ax=ax2,
    label="Normalized overlap"   # <-- questo fa comparire la voce verde
)

offset = 0.322  
for p in ax2.patches:
    p.set_x(p.get_x() + offset)

ax2.set_zorder(ax1.get_zorder() + 1)
ax2.patch.set_visible(False)

ax2.set_ylabel("Normalized overlap", color=color)
ax2.tick_params(axis="y", colors=color)

ax1.set_xlabel("Gene", fontsize=XLABEL_FONTSIZE)
ax1.set_title(args.sample_name + " CNVs summary", fontsize=TITLE_FONTSIZE, fontstyle="italic")

ax1.tick_params(axis="x", rotation=XROT, labelsize=XTICK_FONTSIZE)
for t in ax1.get_xticklabels():
    t.set_fontweight("bold")

ax1.margins(x=XMARGIN)

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()

if ax1.get_legend() is not None:
    ax1.get_legend().remove()
if ax2.get_legend() is not None:
    ax2.get_legend().remove()

leg = ax2.legend(
    h1 + h2,
    l1 + l2,
    fontsize=15,
    frameon=True,
    loc="best"
)

leg.set_zorder(10000)
leg.set_in_layout(False)            
leg.get_frame().set_alpha(0.7)        
leg.get_frame().set_facecolor("white")

fig.savefig(args.summary_folder + "CNV_plot.pdf", bbox_inches="tight")