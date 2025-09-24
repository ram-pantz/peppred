#       Sgourakis Lab
#   Author: Ram Pantula
#   Date: July 1, 2025
#   Email: rpantula@sas.upenn.edu

#!/usr/bin/env python3
import os, sys, subprocess
if len(sys.argv) != 3:
    print("Usage: python genfasta.py <pdb_list_or_file> <output_fasta>")
    sys.exit(1)

inp = sys.argv[1]
out_fa = sys.argv[2]

pdb_dir = "MHC_pdbs"
tmp_dir = "tmp_fastas"
os.makedirs(tmp_dir, exist_ok=True)
os.makedirs(os.path.dirname(out_fa) or ".", exist_ok=True)

# Get list of PDB files
if inp.endswith(".pdb"):
    pdb_files = [inp]
else:
    with open(inp, "r") as f:
        pdb_files = [line.strip() for line in f if line.strip()]

# Start with empty file
open(out_fa, "w").close()

def rewrite_header_and_append(fa_path, new_header, out_file):
    if not os.path.exists(fa_path):
        return False
    with open(fa_path) as f:
        lines = [l.strip() for l in f if l.strip()]
    if not lines:
        return False
    seq = "".join([l for l in lines if not l.startswith(">")])
    with open(out_file, "a") as out:
        out.write(f">{new_header}\n{seq}\n")
    return True

n_ok, n_fail = 0, 0
for pdb in pdb_files:
    fname = os.path.basename(pdb)
    base = fname.replace("_reordered.pdb", "").replace(".pdb", "")
    pdb_path = pdb if os.path.isabs(pdb) else os.path.join(pdb_dir, fname)

    if not os.path.exists(pdb_path):
        print(f"⚠️ Missing PDB: {pdb_path}")
        n_fail += 1
        continue

    faA = os.path.join(tmp_dir, f"{base}_A.fasta")
    faB = os.path.join(tmp_dir, f"{base}_B.fasta")

    # PyMOL script to dump both chains
    pml = f"""
reinitialize
load {pdb_path}, mol
if cmd.count_atoms("mol and chain A") > 0:
    save {faA}, mol and chain A
if cmd.count_atoms("mol and chain B") > 0:
    save {faB}, mol and chain B
quit
"""
    pml_path = os.path.join(tmp_dir, f"{base}.pml")
    with open(pml_path, "w") as f:
        f.write(pml)

    try:
        subprocess.run(["pymol", "-cq", pml_path], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        wroteA = rewrite_header_and_append(faA, f"{base} Chain A", out_fa)
        wroteB = rewrite_header_and_append(faB, f"{base} Chain B", out_fa)
        if wroteA or wroteB:
            n_ok += 1
        else:
            print(f"No chains extracted for {pdb}")
            n_fail += 1
    except Exception as e:
        print(f"PyMOL failed on {pdb}: {e}")
        n_fail += 1

print(f"PDBs OK={n_ok}, Fail={n_fail}. Output → {out_fa}")

