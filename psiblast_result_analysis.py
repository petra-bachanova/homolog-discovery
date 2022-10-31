import pandas as pd
import re


# ........................................FORMAT BLAST.......................................
def merge_and_save(path, df, hs, c, qlen, slen, drop):
    df = pd.merge(left=meta_hs, right=df, left_on="hs_mature_entry", right_on=hs)
    df = pd.merge(left=cono_df, right=df, left_on="c_entry", right_on=c)
    df["scovs"] = (df[qlen] * df["qcovs"].astype(float)) / df[slen]
    if drop:
        df.drop(columns=cols_to_drop).to_csv(f"formatted_data/{path}.csv", index=False)
    else:
        df.to_csv(f"formatted_data/{path}_all.csv", index=False)
    return df


def format_blast(path, hs, c, qlen, slen):
    # path = "220715_psiblast_hsQcDB"
    # hs = "qseqid"
    # c = "sseqid"
    # qlen = "hs_length"
    # slen = "c_length"
    df = pd.read_csv(f"ncbi-blast-2.13.0+/{path}.txt", sep="\t", names=colnames, header=None)
    df = df[df.qcovs > 0]  # qcovs is an integer, tiny qcovs values (for large query protein) are rounded down to 0
    merge_and_save(path, df, hs, c, qlen, slen, drop=True)


def format_psiblast(path, hs, c, qlen, slen):
    path = "221018_psishort_QconoDBhs"
    with open(f"ncbi-blast-2.13.0+/{path}.txt") as f:
        good_lines = []
        lines = [lin for lin in f.readlines() if lin.strip()]
        for lin in lines:
            if lin != "Search has CONVERGED!\n":
                lin = re.sub("\n", "", lin)
                good_lines.append(lin)
    psiblast = pd.DataFrame([lin.split('\t') for lin in good_lines], columns=colnames)
    merge_and_save(path, psiblast, hs, c, qlen, slen, drop=True)


colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']

cols_to_drop = ['qseqid', 'sseqid', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore',
                'c_geneSuperfamily', 'c_cysteineFramewrok', 'c_pharmacologicalFamily', 'c_organismLatin',
                'hs_mature_idx', 'hs_protein names', 'hs_gene Names', 'hs_sequence',
                'hs_post-translational modification', 'hs_involvement in disease']

meta_hs = pd.read_csv("data/221018_metadata_hs.csv")
cono_df = pd.read_csv("data/221017_metadata_c.csv")

format_blast(path="220713_blast_cQhsDB",
             hs="sseqid",
             c="qseqid",
             qlen="c_length",
             slen="hs_length")

format_blast(path="220713_blast_hsQcDB",
             hs="qseqid",
             c="sseqid",
             qlen="hs_length",
             slen="c_length")

format_psiblast(path="220725_pb_cQhsDB_5e3_below30",
                hs="sseqid",
                c="qseqid",
                qlen="c_length",
                slen="hs_length")

format_psiblast(path="220725_pb_hsQcDB_5e3_below30",
                hs="qseqid",
                c="sseqid",
                qlen="hs_length",
                slen="c_length")
