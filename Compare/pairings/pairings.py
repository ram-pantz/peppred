import pandas as pd
from itertools import combinations

df = pd.read_csv("MHC_pdbs/MHCs.csv")[["PDB", "Allele", "Peptide Sequence"]]

same_same, same_diff, diff_diff, all_pairs = [], [], [], []

for (_, a), (_, b) in combinations(df.iterrows(), 2):
    pair = {
        "PDB1": a["PDB"], "PDB2": b["PDB"],
        "Allele1": a["Allele"], "Allele2": b["Allele"],
        "Peptide1": a["Peptide Sequence"], "Peptide2": b["Peptide Sequence"]
    }
    all_pairs.append(pair)

    same_pep = a["Peptide Sequence"] == b["Peptide Sequence"]
    same_hla = a["Allele"] == b["Allele"]

    if same_pep and same_hla:
        same_same.append(pair)
    elif same_pep and not same_hla:
        same_diff.append(pair)
    elif not same_pep and not same_hla:
        diff_diff.append(pair)

# Write all groups
pd.DataFrame(same_same).to_csv("pairings/same_same.csv", index=False)
pd.DataFrame(same_diff).to_csv("pairings/same_diff.csv", index=False)
pd.DataFrame(diff_diff).to_csv("pairings/diff_diff.csv", index=False)
pd.DataFrame(all_pairs).to_csv("pairings/all.csv", index=False)

