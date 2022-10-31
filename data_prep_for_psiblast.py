import re
import numpy as np
import pandas as pd
from datetime import datetime

today = datetime.today().date().strftime("%y%m%d")


def write_seqs_to_fa(f_to_fa, f_name, org):
    """writes entry and sequence of a file to a fasta format"""
    for i in range(len(f_to_fa)):
        with open(f"data/fasta/{today}_{org}_{f_name}.fa", "a") as f:
            f.write(f">{f_to_fa[f'{org}_entry'][i]}\n{f_to_fa[f'{org}_sequence'][i]}\n")


def split_write_seqs(meta, org):
    """splits sequences into those below and above 30aa and uses write_seqs_to_fa fun to save them as fa"""
    below30aa = meta[meta[f'{org}_length'] <= 30].reset_index(drop=True)
    above30aa = meta[meta[f'{org}_length'] > 30].reset_index(drop=True)

    write_seqs_to_fa(f_to_fa=below30aa, f_name='below30', org=org)
    write_seqs_to_fa(f_to_fa=above30aa, f_name='above30', org=org)


def check_repeat_seqs(fa_path):
    """checks how many unique lines are in a fasta files"""
    with open(fa_path, encoding="utf-8") as f:
        lines = f.readlines()
        x = np.array(lines)
        unique_lines = len(np.unique(x))
    print(unique_lines)


# ...................................CONOPEPTIDE DATA...................................
cono_df = pd.read_csv("data/221017_conoserver_structure_protein.csv")
# subset columns and rename
cols_to_subset = ['id', 'name', 'class', 'geneSuperfamily', 'cysteineFramewrok', 'pharmacologicalFamily',
                  'organismLatin', 'sequence', 'sequenceLength', 'UniProt', 'PDB']
new_names = ['c_entry', 'c_name', 'c_class', 'c_geneSuperfamily', 'c_cysteineFramewrok',
             'c_pharmacologicalFamily', 'c_organismLatin', 'c_sequence', 'c_length', 'c_UniProt', 'c_PDB']
cono_df = cono_df[cols_to_subset]
cono_df.columns = new_names
# discard precursor entries
cono_df = cono_df[~cono_df["c_name"].str.contains("precursor") & ~cono_df.duplicated(subset="c_sequence")]
cono_df = cono_df.reset_index(drop=True)

split_write_seqs(meta=cono_df, org='c')
cono_df.to_csv(f"data/{today}_metadata_cono.csv", index=False)


# .....................................HUMAN DATA.......................................
# gff file which contains info about mature chains and peptides of the hs proteins from uniprot
with open("data/uniprot/221017_hs_gff.gff") as file:
    content = file.readlines()
    gff_lines = []
    for line in content:
        if "UniProt" in line:
            gff_lines.append(line)

gff_data = pd.DataFrame([i.split("\t") for i in gff_lines])
gff_data = gff_data.iloc[:, np.r_[0, 2:5, 8]]
gff_data.columns = ["Uniprot_entry", "Seq_attribute", "Start", "Stop", "Info"]
mask = (gff_data['Seq_attribute'] == "Chain") | (gff_data['Seq_attribute'] == "Peptide")
gff_data = gff_data[mask]

mature_entries = []
mature_names = []
for line in gff_data.Info:
    mature_attributes = line.split(";")
    split_key_val = [i.split("=") for i in mature_attributes]
    dict = {item[0]: ''.join(item[1:]) for item in split_key_val}
    mature_entries.append(dict["ID"])
    mature_names.append(dict["Note"])

gff_data["entry"] = mature_entries
gff_data["name"] = mature_names

# same info about hs proteins but in xlsx form with metadata
human_df = pd.read_excel("data/uniprot/221017_hs_excel.xlsx")
human_df.rename(columns={'Entry': 'Uniprot_entry', 'Sequence': 'Full_sequence'}, inplace=True)
human_df = pd.merge(gff_data, human_df, on="Uniprot_entry", how="left")

seq_list = []
len_list = []
idx_list = []
for i in range(len(human_df)):
    seq = human_df['Full_sequence'].loc[i]
    start_idx = int(human_df.Start.loc[i]) - 1
    stop_idx = int(human_df.Stop.loc[i])

    mature_seq = seq[start_idx:stop_idx]
    idx_list.append((start_idx, stop_idx))
    seq_list.append(mature_seq)
    len_list.append(len(mature_seq))

human_df["mature_idx"] = idx_list
human_df["sequence"] = seq_list
human_df["length"] = len_list

human_df = human_df[['Uniprot_entry', 'entry', 'name', 'sequence', 'length']]
human_df.columns = [str('hs_' + str(i[0]).lower() + i[1:]) for i in human_df.columns]
# data.to_csv(f"formatted_data/{today}_metadata_hssecrmat.csv", index=False)
human_df = human_df[human_df.hs_length <= 430]

split_write_seqs(meta=human_df, org='hs')
human_df.to_csv(f"data/{today}_metadata_hs.csv", index=False)
