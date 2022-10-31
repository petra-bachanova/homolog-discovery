import pandas as pd
import json
import re
import os
import shutil
import re
import numpy as np
import glob

# get a list of proteins in folder, e.g. rankings/<name>-ranking_debug
protein_list = [re.sub("-ranking_debug.json", "", i) for i in os.listdir('structure/AlphaFold/output/c/rankings')]
# protein = protein_list[0]

for protein in protein_list:
    # read the debug file, extract the name of the best model
    with open(f"structure/AlphaFold/output/c/rankings/{protein}-ranking_debug.json") as file:
        ranking_df = json.load(file)
        best_model = ranking_df["order"][0]

    # move the good pdb file, e.g. pdbs/<name>-relaxed_<best_model>.pdb to a new folder
    shutil.copy(f"structure/AlphaFold/output/c/pdbs/{protein}-relaxed_{best_model}.pdb",
                f"structure/AlphaFold/output/c/best_models/{protein}-relaxed_{best_model}.pdb")

# .................................
# get a list of all peptide IDs
peptide_ids = glob.glob("structure/AlphaFold/input/c/*.fa", recursive=True)
peptide_ids.extend(glob.glob("structure/AlphaFold/input/hs_1-96/*.fa", recursive=True))
peptide_ids = [re.sub("structure.*?(P.*?).fa", "\\1", i) for i in peptide_ids]

# get a list of structures
with open("structure/AlphaFold/output/best_models/chain_list.txt") as f:
    af_structures = f.read()

af_structures = af_structures.split("\n")
af_structures = [re.sub("-relaxed_.*?$", "", i) for i in af_structures]


def setdiff_sorted(array1, array2, assume_unique=False):
    ans = np.setdiff1d(array1, array2, assume_unique).tolist()
    if assume_unique:
        return sorted(ans)
    return ans


noAFyet = setdiff_sorted(peptide_ids, af_structures)

for peptide in noAFyet:
    shutil.copy(f"structure/AlphaFold/input/all/{peptide}.fa",
                f"structure/AlphaFold/input/noAFyet/{peptide}.fa")


# ---------------------------------------------
# move structures from individual folders
RootDir1 = r'./structure/221025_AF_collab'
TargetFolder = r'./structure/new_pdbs2'

for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
    li = []
    for name in files:
        if '_relaxed_rank_1' in name and name.endswith('.pdb'):
            print("Found relaxed")
            SourceFolder = os.path.join(root, name)
            shutil.copy2(SourceFolder, TargetFolder)
            li.append('relaxed')
    if len(li) == 0:
        for name in files:
            if '_unrelaxed_rank_1' in name and name.endswith('.pdb'):
                print('unrelaxed')
                SourceFolder = os.path.join(root, name)
                shutil.copy2(SourceFolder, TargetFolder)

