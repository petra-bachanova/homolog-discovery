import matplotlib.ticker
from datetime import datetime
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def plot_heatmap(dat, dim, name, cluster, dendro):
    pivot = pd.pivot_table(dat, index="hs_name", columns="c_name", values="TM1")
    pivot = pivot.fillna(0)
    if not cluster:
        kws = dict(cbar_kws=cbar_dictionary, figsize=dim, col_cluster=False, row_cluster=False)
    else:
        kws = dict(cbar_kws=cbar_dictionary, figsize=dim)
    g = sns.clustermap(pivot, cmap=cmap, norm=norm, xticklabels=True, yticklabels=True, **kws)
    g.ax_cbar.set_title('TM1 score')
    ax = g.ax_heatmap
    ax.set_ylabel("")
    ax.set_xlabel("")
    if not dendro:
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
    # g.ax_cbar.set_position([0.3, 0.92, 0.6, 0.3])  # settings for colorbar
    g.ax_cbar.set_position([0.02, 0.9, 0.12, 0.02])  # settings for heatmap with all
    # g.ax_cbar.set_position([0.02, 0.95, 0.12, 0.02])
    plt.savefig(f"graphs/Fig3_heatmaps/{name}.png", dpi=300)


def plot_TM(dat, name, quadrants, marker_s, colour_by, cbar_title):
    ax = sns.scatterplot(x=dat['TM2'], y=dat['TM1'], hue=dat[colour_by], s=marker_s, palette='Spectral_r', linewidth=0)
    norm = plt.Normalize(dat['lenDiff'].min(), dat[colour_by].max())
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
    if quadrants:
        ax.axhline(0.4, c='#616075', linewidth=1, linestyle='--')
        ax.axvline(0.4, c='#616075', linewidth=1, linestyle='--')
    sm.set_array([])
    ax.get_legend().remove()
    cbar = ax.figure.colorbar(sm)
    cbar.set_label(cbar_title)
    plt.savefig(f"graphs/Fig2_psi_tmalign/{today}_TM_{name}.png", dpi=300)


def plot_TM_psi_metrics(metric, cbar_title, palette):
    # metric = 'Qcono_pident'
    mask = tm_psi[metric].isna()
    plt.figure(figsize=(5, 4), dpi=100)
    sns.scatterplot(x=tm_psi[mask]['TM2'], y=tm_psi[mask]['TM1'], color='#CAB0B5', s=15, linewidth=0)
    ax = sns.scatterplot(x=tm_psi['TM2'], y=tm_psi['TM1'], hue=tm_psi[metric], s=15, palette=palette,
                         linewidth=0)
    norm = plt.Normalize(tm_psi[metric].min(), tm_psi[metric].max())
    sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
    ax.axhline(0.4, c='#616075', linewidth=1, linestyle='--')
    ax.axvline(0.4, c='#616075', linewidth=1, linestyle='--')
    sm.set_array([])
    ax.get_legend().remove()
    cbar = ax.figure.colorbar(sm)
    cbar.set_label(f'{cbar_title} as query search')
    plt.savefig(f"graphs/Fig2_psi_tmalign/{today}_TM_psi_{metric}.png", dpi=300)
    plt.close()


today = datetime.now().strftime('%y%m%d')
# imports
meta_hs = pd.read_csv("./data/221018_metadata_hs.csv")
meta_c = pd.read_csv("./data/221017_metadata_cono.csv")
tm_data = pd.read_csv("./structure/221025_final pdbs/221025_tmalign.txt", sep="\t")
tm_data.rename(columns={"#PDBchain1": "PDBchain1"}, inplace=True)

# strip chain names (now file paths) to protein entries
for col in ["PDBchain1", "PDBchain2"]:
    tm_data[col] = tm_data[col].str.replace(r'.[/]', '', regex=True)
    tm_data[col] = tm_data[col].str.replace(r'(P[0-9]{5}).*?$', '\\1', regex=True)
    tm_data[col] = tm_data[col].str.replace(r'(PRO_[0-9]{10}).*?$', '\\1', regex=True)

# get alignment comparisons between conopeptides and human proteins rather than within either group
mask1 = (tm_data["PDBchain1"].str.contains("^P[0-9]")) & (tm_data["PDBchain2"].str.contains("^PRO"))
tm_data = tm_data[mask1]
tm_data.rename(columns={'PDBchain1': 'cono', 'PDBchain2': 'hs'}, inplace=True)
cols_c = ['c_entry', 'c_name', 'c_length', 'c_sequence']
cols_hs = ['hs_entry', 'hs_name', 'hs_length', 'hs_sequence']
tm_data = pd.merge(tm_data, meta_c[cols_c], left_on="cono", right_on="c_entry")
tm_data = pd.merge(tm_data, meta_hs[cols_hs], left_on="hs", right_on="hs_entry")
# tm_data.to_csv('./data/221025_tm_data_no_threshold.csv', index=False)
# plot TM1 vs TM2, colour by difference in lengths
tm_data['lenDiff'] = abs(tm_data['c_length'] - tm_data['hs_length'])
plot_TM(dat=tm_data, name='lenDiff', quadrants=True, marker_s=8,
        colour_by='lenDiff', cbar_title='Absolute difference in homolog lengths')
plot_TM(dat=tm_data, name='RMSD', quadrants=True, marker_s=8, colour_by='RMSD', cbar_title='RMSD')

psi = pd.read_csv('data/221028_psiblast_metrics.csv')
tm_psi = pd.merge(tm_data, psi, on=['cono', 'hs'], how='outer')

metric_li = ['Qcono_pident', 'Qcono_evalue', 'Qcono_qcovs', 'Qhs_pident', 'Qhs_evalue', 'Qhs_qcovs']
title_li = ['% identity for Cono', 'Eval for Cono', 'Query coverage for Cono',
            '% identity for Hs protein', 'Eval for Hs protein', 'Query coverage for Hs protein']
pal = ['RdYlGn_r'] * 3 + ['RdYlBu_r'] * 3

for i in range(6):
    plot_TM_psi_metrics(metric=metric_li[i], cbar_title=title_li[i], palette=pal[i])

color_palette = sns.diverging_palette(250, 0, as_cmap=True)

# TM scores where conopeptide aligned to human > 0.5 and human aligned to cono is > 0.5
mask2 = (tm_data["TM1"] > 0.4) & (tm_data["TM2"] > 0.4)
good_data = tm_data[mask2].copy()
long_name = 'Platelet-derived growth factor D%2C receptor-binding form'
short_name = 'Platelet-derived GF-D receptor-binding form'
good_data.loc[good_data.hs_name == long_name, 'hs_name'] = short_name
good_data.to_csv('./data/221025_tm_data_tm1-04_tm2-04.csv', index=False)

# plot TM1 vs TM2, colour by difference in lengths
plot_TM_lenDiff(dat=good_data, name='good_data', quadrants=False, marker_s=15)

# tropic hormones and growth factors
tropic_li = ['Glycoprotein hormone beta-5 prepropeptide',
             'Glycoprotein hormone alpha-2 prepropeptide',
             'thyrostimulin-beta 5', 'thyrostimulin-alpha 2']
mask_tropic = good_data["c_name"].isin(tropic_li)
tropic = good_data[mask_tropic]
# conkunitzins
mask_conk = good_data["c_name"].str.contains("Conkunitz")
conk = good_data[mask_conk]
# insulins
mask_ins = good_data["hs_name"].str.contains("Ins")
ins = good_data[mask_ins]
# smaller groups
rest = good_data[~(mask_conk | mask_tropic | mask_ins)]

# heatmap parameters
cmap = sns.color_palette("viridis", as_cmap=True)
norm = mcolors.LogNorm()
sns.set(font_scale=0.8)
cbar_dictionary = dict(ticks=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                       format=matplotlib.ticker.ScalarFormatter(useMathText=False),
                       orientation='horizontal')

plot_heatmap(dat=good_data, dim=(18, 10), name=f'{today}_heatmap_all_supplementary.png', cluster=True, dendro=True)
plot_heatmap(dat=tropic, dim=(4, 6.66), name=f'{today}_tropic.png', cluster=True, dendro=False)
plot_heatmap(dat=conk, dim=(10.57, 2.63), name=f'{today}_conk.png', cluster=True, dendro=False)
plot_heatmap(dat=ins, dim=(3, 2), name=f'{today}_ins.png', cluster=False, dendro=False)
plot_heatmap(dat=rest, dim=(11.56, 7), name=f'{today}_rest.png', cluster=True, dendro=False)
plot_heatmap(dat=ins, dim=(3, 4), name=f'{today}_colorbar.png', cluster=False, dendro=False)
