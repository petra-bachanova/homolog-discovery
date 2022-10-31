from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import glob
import re


def format_psiblast_output(psi1, psi2, query, subject):
    interesting_cols = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs']
    psi = pd.concat([psiblast_frames[psi1], psiblast_frames[psi2]])[interesting_cols]
    psi = psi.drop_duplicates(subset=['qseqid', 'sseqid'], keep='last')
    psi.columns = [re.sub('^', f'Q{query}_', i) for i in psi.columns]
    psi.rename(columns={f'Q{query}_qseqid': f'{query}', f'Q{query}_sseqid': f'{subject}'}, inplace=True)
    return psi


today = datetime.now().strftime('%y%m%d')
psiblast_files = glob.glob('./ncbi-blast-2.13.0+/*psi*.txt')
colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']

all_proteins = []
psiblast_frames = []
for file in psiblast_files:
    # file = files[0]
    with open(file) as f:
        good_lines = []
        lines = [lin for lin in f.readlines() if lin.strip()]
        for lin in lines:
            if lin != "Search has CONVERGED!\n":
                lin = re.sub("\n", "", lin)
                good_lines.append(lin)
    psiblast_frame = pd.DataFrame([lin.split('\t') for lin in good_lines], columns=colnames)
    psiblast_frames.append(psiblast_frame)
    li = psiblast_frame['qseqid'].tolist()
    li.extend(psiblast_frame['sseqid'].tolist())
    all_proteins.extend(li)

# psiblast_frames: psishort_QconoDBhs.txt, psishort_QhsDBcono.txt, psi_QconoDBhs.txt, psi_QhsDBcono.txt'
QconoDBhs = format_psiblast_output(psi1=0, psi2=2, query='cono', subject='hs')
QhsDBcono = format_psiblast_output(psi1=1, psi2=3, query='hs', subject='cono')

merged = pd.merge(QconoDBhs, QhsDBcono, on=['hs', 'cono'], how='outer')
# merged.to_csv(f'data/{today}_psiblast_metrics.csv', index=False)

# len(merged.dropna()) is the length of the shared elements (otherwise one of the metrics for Qhs/Qcono would be na)
plt.rcParams.update({'font.size': 15})
plt.figure(figsize=(5, 4))
venn2(subsets=(len(QconoDBhs), len(QhsDBcono), len(merged.dropna())),
      set_labels=('Cono', ''),
      set_colors=('#46ac5c', '#699fc9'), alpha=0.7)
plt.savefig(f'graphs/Fig2_psi_tmalign/{today}_homolog_pairs_venn.png', dpi=300)


all_proteins = set(all_proteins)
human = [i for i in all_proteins if 'PRO' in i]
len(human)  # 181
cono = [i for i in all_proteins if 'PRO' not in i]
len(cono)  # 362

# what human hits have been found?
human_meta = pd.read_csv('data/221018_metadata_hs.csv')
human_hits_meta = human_meta[human_meta['hs_entry'].isin(human)]
human_hits_meta = human_hits_meta.sort_values(by='hs_name')

#
#
# get a list of proteins the structures of which have already been resolved
prot_with_structure = glob.glob(r'./structure/*')
prot_with_structure = [i.replace("./structure\\", "") for i in prot_with_structure]
prot_with_structure = [re.sub('(P.*?)-relaxed.*?$', '\\1', i) for i in prot_with_structure]  # for human peptides
prot_with_structure = [re.sub('(P[0-9]{5})_.*?$', '\\1', i) for i in prot_with_structure]  # for conopeptides

homologs_with_structure = set(all_proteins).intersection(set(prot_with_structure))  # 309
len(homologs_with_structure)
to_run = [i for i in all_proteins if i not in homologs_with_structure]
to_run.sort()

human_to_run = human_meta[human_meta['hs_entry'].isin(to_run)]
human_to_run[['hs_entry', 'hs_name', 'hs_sequence']].to_csv('data/221023_human_to_run.csv', index=False)

cono_meta = pd.read_csv('data/221017_metadata_cono.csv')
cono_to_run = cono_meta[cono_meta['c_entry'].isin(to_run)]
cono_to_run[['c_entry', 'c_name', 'c_sequence']].to_csv('data/221023_cono_to_run.csv', index=False)
