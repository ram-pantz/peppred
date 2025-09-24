#       Sgourakis Lab
#   Author: Ram Pantula
#   Date: July 1, 2025
#   Email: rpantula@sas.upenn.edu

import os
import glob
import numpy as np
import pandas as pd
import joblib
import itertools
import matplotlib.pyplot as plt

from misc.common import *
from Bio.PDB import PDBParser, NeighborSearch, is_aa  # <-- required

# CONFIGURATION
PAIRS_CSV       = f"{dbpath}/pairings/same_diff.csv"
SCORING_DIR     = f"{dbpath}/scorings"
PDB_DIR         = f"{dbpath}/MHC_pdbs"
ANCHOR_CSV      = f"{dbpath}/anchor_class/anchor_class.csv"
DIHEDRAL_STEM   = f"{dbpath}/peptide_dihedrals/peptide_dihedrals_"
SIM_THRESHOLD   = 1.5

# Dihedral column names
dihed_cols = ["phi_4","psi_4","phi_5","psi_5","phi_6","psi_6","phi_7","psi_7"]


# --------------------------
# Feature builder
# --------------------------
def read_central_dihedral(file_stem: str, pdbs_list, anchor_class: int):
    """
    Load dihedrals for a given anchor class from ONE CSV:
    peptide_dihedrals_{anchor_class+2}.csv
    Returns dict: {pdbid: np.array([... 8 angles in radians ...])}
    """
    dihed_dict = {}
    file_path = f"{file_stem}{anchor_class + 2}.csv"

    if not os.path.exists(file_path):
        print(f"DEBUG: Dihedral CSV file not found: {file_path}")
        return dihed_dict

    print(f"DEBUG: Found dihedral file: {file_path}")
    try:
        df = pd.read_csv(file_path)
        df = df[df["pdbid"].isin(pdbs_list)].copy()

        if any(col not in df.columns for col in dihed_cols) or df.empty:
            print(f"DEBUG: Required columns missing or no matching PDBs in {file_path}. Skipping.")
            return dihed_dict

        for _, row in df.iterrows():
            dihedrals = row[dihed_cols].to_numpy(dtype=float, copy=False).ravel()
            dihed_dict[str(row["pdbid"])] = dihedrals

        print(f"DEBUG: Successfully loaded dihedral data for {len(dihed_dict)} PDBs.")
    except Exception as e:
        print(f"DEBUG: Error loading {file_path}: {e}. Skipping.")
    return dihed_dict


def circ_abs_diff_rad(a, b):
    """Minimal absolute angular difference for radians: returns value in [0, pi]."""
    x = (a - b + np.pi) % (2*np.pi) - np.pi
    return np.abs(x)


def _find_pdb_path(pdbid):
    """Try common filename patterns in PDB_DIR."""
    candidates = [
        os.path.join(PDB_DIR, f"{pdbid}.pdb"),
        os.path.join(PDB_DIR, f"{pdbid}_base.pdb"),
        os.path.join(PDB_DIR, f"{pdbid}_reordered.pdb"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    g = glob.glob(os.path.join(PDB_DIR, f"{pdbid}*.pdb"))
    return g[0] if g else None


def _chainA_min_distances_to_chainB(pdb_path):
    """
    For a PDB with chains 'A' (HLA) and 'B' (peptide),
    return a list of minimal heavy-atom distances (Å) from each AA in chain A to any atom in chain B.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("x", pdb_path)
        model = structure[0]
        chain_a = model["A"]
        chain_b = model["B"]
    except Exception:
        return []

    atoms_b = [atom for res in chain_b for atom in res if getattr(atom, "element", "") != "H"]
    if not atoms_b:
        return []

    ns = NeighborSearch(list(atoms_b))
    dists = []
    for res in chain_a:
        if not is_aa(res):
            continue
        min_d = float("inf")
        for atom in res:
            if getattr(atom, "element", "") == "H":
                continue
            neighbors = ns.search(atom.coord, 12.0, level='A')  # 'A' = Atom level
            if neighbors:
                for nb in neighbors:
                    d = atom - nb
                    if d < min_d:
                        min_d = d
        dists.append(min_d if np.isfinite(min_d) else 1e6)  # large if none found
    return dists


def _bin_weight_from_distance(d):
    """Piecewise weights: closer to chain B → higher weight."""
    if d <= 3.0:   return 1.0
    if d <= 6.0:   return 0.7
    if d <= 9.0:   return 0.4
    return 0.2


def _compute_hla_weights(pdb_path, hla_len):
    """
    Returns a weight vector (length = hla_len) for chain A residues, based on distance to chain B.
    If pdb missing or length mismatch, fall back gracefully to a low, uniform weight.
    """
    if (pdb_path is None) or (not os.path.exists(pdb_path)):
        return np.full(hla_len, 0.4, dtype=float)

    dists = _chainA_min_distances_to_chainB(pdb_path)
    if not dists:
        return np.full(hla_len, 0.4, dtype=float)

    w = np.array([_bin_weight_from_distance(d) for d in dists], dtype=float)

    # Align to expected length (hla_len). If different, truncate or pad with last value.
    if len(w) == hla_len:
        return w
    out = np.full(hla_len, w[-1] if len(w) else 0.4, dtype=float)
    n = min(hla_len, len(w))
    out[:n] = w[:n]
    return out


def extract_features_for_prediction(df_pairs: pd.DataFrame, dihed_dict: dict):
    """
    Returns:
        X_base : [n, 16] (circ-deltas + pLDDT summaries)
    """
    X = []
    hla_slice = slice(0, 180)   # HLA portion of pLDDT
    pep_slice = slice(180, 189) # last 9 residues = peptide

    # simple counters for visibility
    dropped = {"no_dihed":0, "npz_missing":0, "plddt_short":0}

    for _, r in df_pairs.iterrows():
        PDB1, PDB2 = str(r.PDB1), str(r.PDB2)
        if PDB1 not in dihed_dict or PDB2 not in dihed_dict:
            dropped["no_dihed"] += 1
            continue

        # circular dihedral deltas (radians)
        dA = np.asarray(dihed_dict[PDB1], dtype=float)
        dB = np.asarray(dihed_dict[PDB2], dtype=float)
        if dA.shape != dB.shape or dA.size != len(dihed_cols):
            dropped["no_dihed"] += 1
            continue
        circ_delta = circ_abs_diff_rad(dA, dB)

        # pLDDT vectors (strict key names per archive)
        try:
            with np.load(os.path.join(SCORING_DIR, f"{PDB1}.npz"), allow_pickle=True) as npzA:
                plddt_A = npzA[f"outfile_{PDB1}_model_1_model_2_ptm_ft_plddt.npy"]
            with np.load(os.path.join(SCORING_DIR, f"{PDB2}.npz"), allow_pickle=True) as npzB:
                plddt_B = npzB[f"outfile_{PDB2}_model_1_model_2_ptm_ft_plddt.npy"]
        except (KeyError, FileNotFoundError):
            dropped["npz_missing"] += 1
            continue

        plddt_A = np.asarray(plddt_A, dtype=float)
        plddt_B = np.asarray(plddt_B, dtype=float)
        if plddt_A.ndim != 1 or plddt_B.ndim != 1 or len(plddt_A) < 189 or len(plddt_B) < 189:
            dropped["plddt_short"] += 1
            continue

        # distance-weighted HLA pLDDT vectors
        wA = _compute_hla_weights(_find_pdb_path(PDB1), len(plddt_A[hla_slice]))
        wB = _compute_hla_weights(_find_pdb_path(PDB2), len(plddt_B[hla_slice]))
        wA = wA / (wA.sum() + 1e-8)
        wB = wB / (wB.sum() + 1e-8)
        hla_wavg_A = float(np.dot(plddt_A[hla_slice], wA))
        hla_wavg_B = float(np.dot(plddt_B[hla_slice], wB))

        # peptide (last 9) + unweighted HLA means (context)
        pep_avg_A = float(np.mean(plddt_A[pep_slice]))
        pep_avg_B = float(np.mean(plddt_B[pep_slice]))
        hla_avg_A = float(np.mean(plddt_A[hla_slice]))
        hla_avg_B = float(np.mean(plddt_B[hla_slice]))

        # pairwise combos
        hla_wavg_absdiff = abs(hla_wavg_A - hla_wavg_B)
        hla_wavg_mean    = 0.5 * (hla_wavg_A + hla_wavg_B)
        pep_avg_absdiff  = abs(pep_avg_A - pep_avg_B)
        pep_avg_mean     = 0.5 * (pep_avg_A + pep_avg_B)

        X.append(np.hstack([
            circ_delta,
            [hla_wavg_absdiff, hla_wavg_mean, pep_avg_absdiff, pep_avg_mean,
             hla_avg_A, hla_avg_B, pep_avg_A, pep_avg_B]
        ]))

    if dropped["no_dihed"] or dropped["npz_missing"] or dropped["plddt_short"]:
        print(f"DEBUG: Dropped pairs -> no_dihed={dropped['no_dihed']}, npz_missing={dropped['npz_missing']}, plddt_short={dropped['plddt_short']}")

    if not X:
        return np.empty((0, 16))

    return np.vstack(X).astype(float)

def produce_final_outputs(
    df_with_preds: pd.DataFrame,
    outfile: str = "final_outputs.txt",
    peptide_cols: tuple | list | None = None,   # e.g. ("Pep1","Pep2") or ("peptide",)
    hla_cols: tuple | list | None = None        # e.g. ("HLA1","HLA2")
):
    """
    Create a text file:

    PEPTIDE
    HLA1 - HLA2
    HLA1 - HLA2

    • Only rows with bucket == "Similar".
    • If two peptide columns exist (Pep1/Pep2), the peptide is taken from both sides.
    • You may override auto-detection via peptide_cols / hla_cols.
    """

    if "bucket" not in df_with_preds.columns:
        raise ValueError("DataFrame must include 'bucket'.")

    df = df_with_preds[df_with_preds["bucket"] == "Similar"].copy()
    if df.empty:
        with open(outfile, "w") as f:
            f.write("# No Similar pairs found.\n")
        print(f"Wrote empty output (no Similar pairs) to: {outfile}")
        return

    cols_lower = {c: c.lower() for c in df.columns}

    # ---- Detect peptide columns ----
    if peptide_cols is None:
        # candidates by regex-ish heuristics
        pep_like = [c for c in df.columns
                    if ("pep" in c.lower() or "seq" in c.lower() or "peptid" in c.lower())]
        # prefer exactly-two columns (Pep1/Pep2 style)
        if len(pep_like) >= 2:
            # try to pick the two with trailing 1/2, else first two alphabetically
            p1s = [c for c in pep_like if c.lower().endswith(("1", "_1"))]
            p2s = [c for c in pep_like if c.lower().endswith(("2", "_2"))]
            if p1s and p2s:
                peptide_cols = (p1s[0], p2s[0])
            else:
                peptide_cols = tuple(sorted(pep_like)[:2])
        elif len(pep_like) == 1:
            peptide_cols = (pep_like[0],)
        else:
            # last chance: sometimes peptide sequence is embedded in PDB-derived names; bail gracefully
            raise ValueError(
                f"Could not auto-detect peptide columns. Available columns: {list(df.columns)}"
            )
    else:
        # ensure tuple-ish
        peptide_cols = tuple(peptide_cols)

    # ---- Detect HLA columns ----
    if hla_cols is None:
        hla_like = [c for c in df.columns if ("hla" in c.lower() or "allele" in c.lower())]
        if len(hla_like) >= 2:
            # prefer columns with 1/2 suffix; else first two sorted for determinism
            h1s = [c for c in hla_like if c.lower().endswith(("1", "_1"))]
            h2s = [c for c in hla_like if c.lower().endswith(("2", "_2"))]
            if h1s and h2s:
                hla_cols = (h1s[0], h2s[0])
            else:
                hla_cols = tuple(sorted(hla_like)[:2])
        else:
            raise ValueError(
                f"Could not auto-detect two HLA/allele columns. Available columns: {list(df.columns)}"
            )
    else:
        hla_cols = tuple(hla_cols)
        if len(hla_cols) != 2:
            raise ValueError("hla_cols must be a 2-tuple/list like ('HLA1','HLA2').")

    # ---- Normalize + build long view by peptide ----
    key_cols = list(set(peptide_cols) | set(hla_cols))
    df = df.dropna(subset=key_cols).copy()
    for c in key_cols:
        df[c] = df[c].astype(str)

    # Build (peptide, pair) rows:
    rows = []
    if len(peptide_cols) == 1:
        pcol = peptide_cols[0]
        for _, r in df.iterrows():
            rows.append((r[pcol], f"{r[hla_cols[0]]} - {r[hla_cols[1]]}"))
    else:
        p1, p2 = peptide_cols[:2]
        for _, r in df.iterrows():
            pair_str = f"{r[hla_cols[0]]} - {r[hla_cols[1]]}"
            rows.append((r[p1], pair_str))
            rows.append((r[p2], pair_str))

    # ---- Write grouped output ----
    from collections import defaultdict
    by_pep = defaultdict(set)
    for pep, pair in rows:
        by_pep[pep].add(pair)

    with open(outfile, "w") as f:
        for pep in sorted(by_pep):
            f.write(f"{pep}\n")
            for pair in sorted(by_pep[pep]):
                f.write(f"{pair}\n")
            f.write("\n")

    print(f"Wrote final outputs to: {outfile}")


# ===========================
# REVISED MAIN FUNCTION FOR PREDICTION
# ===========================
def main():
    # Load your new pairs
    df_new = pd.read_csv(PAIRS_CSV)

    # Load the pre-trained model and scaler (only here; no duplicate globals)
    model = joblib.load("model.pkl")
    scaler = joblib.load("scaler.pkl")

    # Load the thresholds
    thr_low = float(np.load("thr_low.npy").reshape(-1)[0])
    thr_high = float(np.load("thr_high.npy").reshape(-1)[0])

    # Get a list of all unique PDBs
    pdbs = list(set(df_new["PDB1"].astype(str)).union(set(df_new["PDB2"].astype(str))))

    # Load the dihedral data (anchor_class=7 → file _9.csv)
    dihed_dict = read_central_dihedral(DIHEDRAL_STEM, pdbs, 7)

    if not dihed_dict:
        print("Error: Dihedral data not found or could not be loaded. Exiting.")
        return

    # Keep only pairs that have dihedrals
    df_filtered = df_new[df_new.PDB1.astype(str).isin(dihed_dict) & df_new.PDB2.astype(str).isin(dihed_dict)].reset_index(drop=True)

    # Extract features
    X_new = extract_features_for_prediction(df_filtered, dihed_dict)

    if X_new.size == 0:
        print("No usable pairs found for prediction.")
        return

    # Scale + predict
    X_new_scaled = scaler.transform(X_new)
    y_proba = model.predict_proba(X_new_scaled)[:, 1]

    # Buckets via thresholds
    predictions = np.where(
        y_proba >= thr_high, "Similar",
        np.where(y_proba >= thr_low, "Uncertain", "Dissimilar")
    )

    # Save results
    df_out = df_filtered.copy()
    df_out["confidence"] = y_proba
    df_out["bucket"] = predictions

    print(df_out.head())
    df_out.to_csv("predictions.csv", index=False)
    print("Predictions saved to predictions.csv")
    
    produce_final_outputs(
    df_out,
    outfile="final_outputs.txt",
    peptide_cols=("Peptide1","Peptide2"),
    hla_cols=("Allele1","Allele2")
    )


if __name__ == "__main__":
    main()

