import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from datetime import datetime
from sklearn.preprocessing import MinMaxScaler


def numbers_out(df):
    return (f"There are {len(df['c_entry'].unique())} unique conopeptides, "
            f"{len(df['hs_mature_entry'].unique())} unique human peptides "
            f"and {df.groupby(['c_entry', 'hs_mature_entry']).ngroups} unique combinations.")


def format_blast_output(df, entry, metadata, tuple_of):
    output = df[entry].value_counts()
    output = pd.DataFrame({entry: output.index, 'count': output.values})
    output = pd.merge(left=metadata, right=output, on=entry)
    homologs = [tuple(df[df[entry] == e][tuple_of].unique()) for e in output[entry]]
    output["homologue_list"] = homologs
    return output


def validate_known_homologs(df, groupby, top_n):
    cono_list = ["contulakin", "conorfamide", "conoNPY", "conopressin", "conoCAP",
                 "con-ins", "thyrostimulin", "conorphin", "MiXXVIIA"]
    hs_list = ["neurotensin", "pro-FMRFamide-related", "neuropeptide", "vasopressin",
               "insulin", "thyrotropin", "dynorphin", "granulin A", "oxytocin"]
    df = df[df['hs_mature_name'].str.contains('|'.join(hs_list), na=False, flags=re.IGNORECASE)]
    df = df[df['c_name'].str.contains('|'.join(cono_list), na=False, flags=re.IGNORECASE)]

    metrics = ['pident', 'qcovs', 'evalue', 'scovs']
    df = df.sort_values(by=metrics, ascending=[False, False, True, False]).groupby(groupby).head(top_n)
    return df


def homolog_heatmap(df, name, colour_by, fillna):
    # df = test_hsdb
    pivot = pd.pivot_table(df, index="c_name", columns="hs_mature_name", values=colour_by)
    pivot = pivot.fillna(fillna)

    cmap = plt.cm.viridis
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[-1] = (1.0, 1.0, 1.0, 1.0)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('mcm', cmaplist, cmap.N)

    sns.set(font_scale=0.5)
    fontsize_pt = 7
    dpi = 72.27
    norm = mcolors.LogNorm()
    matrix_height_pt = fontsize_pt * pivot.shape[0]
    matrix_height_in = matrix_height_pt / dpi
    top_margin = 0.04  # in percentage of the figure height
    bottom_margin = 0.04  # in percentage of the figure height
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
    sns.clustermap(pivot, cmap=cmap, figsize=(5, figure_height), norm=norm, xticklabels=True, yticklabels=True)
    # plt.savefig(f"results/{name}.png", dpi=600)


today = datetime.today().date().strftime("%y%m%d")
meta_hs = pd.read_csv("formatted_data/220713_meta_hs.csv")
meta_c = pd.read_csv("formatted_data/220713_meta_c.csv")

# ..................................PREP FOR ALPHA FOLD..................................
hsdb_below30 = pd.read_csv("formatted_data/220725_pb_cQhsDB_5e3_below30.csv")
hsdb_above30 = pd.read_csv("formatted_data/220725_pb_cQhsDB_5e3_above30.csv")

cdb_below30 = pd.read_csv("formatted_data/220725_pb_hsQcDB_5e3_below30.csv")
cdb_above30 = pd.read_csv("formatted_data/220725_pb_hsQcDB_5e3_above30.csv")

df_list = [hsdb_below30, hsdb_above30, cdb_below30, cdb_above30]

for i in df_list:
    print(numbers_out(i))

hs = []
c = []
for i in df_list:
    hs.extend(i.hs_mature_name)
    c.extend(i.c_name)
hs = np.unique(np.array(hs))
c = np.unique(np.array(c))

plt.hist(meta_hs[meta_hs["hs_mature_entry"].isin(hs)]["hs_length"], bins=30, color="#515", alpha=0.7)
plt.hist(meta_c[meta_c["c_entry"].isin(c)]["c_length"], bins=30, color="#0B5E55", alpha=0.7)
plt.savefig("./results/Distribution of hit lengths", dpi=300)


# make fasta files of interesting hs peptides
def write_AlphaFold_files(meta, peptide_list, entry, sequence, length, list_name):
    # list of peptide names sorted by ascending length (will be run on AF in this order)
    meta_in_list = meta[meta[entry].isin(peptide_list)]
    meta_in_list = meta_in_list.sort_values(by=length)
    with open(f"./structure/AlphaFold/{list_name}.txt", "w") as file:
        for index in range(len(meta_in_list)):
            file.write(f"'{meta_in_list[entry].values[index]}', ")

    for peptide in peptide_list:
        seq = meta.loc[meta[entry] == peptide][sequence].values[0]
        with open(f"./structure/AlphaFold/input/c/{peptide}.fa", "a") as file:
            file.write(f">{peptide}\n{seq}\n")


write_AlphaFold_files(meta=meta_c,
                      peptide_list=c,
                      entry="c_entry",
                      sequence="c_sequence",
                      length="c_length",
                      list_name="conopeptide_list")

# ..................................EXPLORING RELATIONSHIPS BETWEEN BLAST METRICS..................................
data_cdb = pd.DataFrame()
data_hsdb = pd.DataFrame()

cols = ['pident', 'qcovs', 'evalue', 'scovs', 'hs_length', 'c_length']

g = sns.pairplot(data_cdb[cols])
plt.subplots_adjust(bottom=0.05)
g.fig.suptitle("Pairplot for hs queries of cono db", y=1.001)
plt.close(g)
plt.savefig("Pairplot_psiblast_hsQcDB.jpeg", dpi=300)

g = sns.pairplot(data_hsdb[cols])
plt.subplots_adjust(bottom=0.05)
g.fig.suptitle("Pairplot for cono queries of hs db", y=1.001)
plt.savefig("Pairplot_psiblast_cQhsDB.jpeg", dpi=300)

scaler = MinMaxScaler()
scaled = scaler.fit_transform(data_cdb[cols])
sns.heatmap(pd.DataFrame(scaled, columns=cols).corr(), annot=True)
plt.savefig("Correlation_psi_hsQcDB.jpeg", dpi=300)

scaled = scaler.fit_transform(data_hsdb[cols])
sns.heatmap(pd.DataFrame(scaled, columns=cols).corr(), annot=True)
plt.savefig("Correlation_psi_cQhsDB.jpeg", dpi=300)

# ..................................VALIDATING THE SEARCH WITH KNOWN HOMOLOGS..................................
known_hsdb = validate_known_homologs(df=data_hsdb, groupby='hs_mature_name', top_n=5)
known_cdb = validate_known_homologs(df=data_cdb, groupby='hs_mature_name', top_n=5)

writer = pd.ExcelWriter('known_homologs.xlsx', engine='xlsxwriter')
known_hsdb.to_excel(writer, sheet_name='known_hsDB_conoQ')
known_cdb.to_excel(writer, sheet_name='known_conoDB_hsQ')


# ..................................................DISCOVERY..................................................
def make_discovery(df, scovs, qcovs, pident, head, group):
    df = df[
        (df.scovs > scovs) & (df.qcovs > qcovs) & (df.pident > pident)]
    df = df.sort_values(by=['hs_mature_name', 'pident', 'qcovs', 'evalue', 'scovs'],
                        ascending=[False, False, False, True, False]).groupby(group).head(head)
    print(numbers_out(df))
    return df


# test_hsdb = make_discovery(df=data_hsdb, scovs=50, qcovs=20, pident=10, head=100, group="c_name")
# test_cdb = make_discovery(df=data_cdb, scovs=20, qcovs=50, pident=10, head=100, group="hs_mature_name")

test_hsdb = data_hsdb[(data_hsdb.scovs > 50) & (data_hsdb.qcovs > 20) & (data_hsdb.pident > 10)]
test_cdb = data_cdb[(data_cdb.scovs > 50) & (data_cdb.qcovs > 20) & (data_cdb.pident > 10)]

homolog_heatmap(df=test_cdb, name="cdb_heatmap", colour_by="evalue", fillna=1)
homolog_heatmap(df=test_hsdb, name="hsdb_heatmap", colour_by="evalue", fillna=1)
