#!/usr/bin/env python3
import os
import sys
import shutil

def main():
    if len(sys.argv) < 2:
        print("Usage: python move_folders.py <root_directory>")
        sys.exit(1)

    root = os.path.abspath(sys.argv[1])

    afft_dir = os.path.join(root, "AFFT-HLA3DB")
    input_seq_dir = os.path.join(afft_dir, "input_seq")
    outfiles_dir = os.path.join(afft_dir, "outfiles")
    mhc_pdbs_dir = os.path.join(afft_dir, "MHC_pdbs")

    # Ensure needed directories exist
    os.makedirs(outfiles_dir, exist_ok=True)
    os.makedirs(mhc_pdbs_dir, exist_ok=True)

    # Step 1: Move matching folders into outfiles
    for fname in os.listdir(input_seq_dir):
        if not fname.endswith("_seq.txt"):
            continue

        targetfolder = fname.replace("_seq.txt", "")
        src = os.path.join(afft_dir, targetfolder)
        dst = os.path.join(outfiles_dir, targetfolder)

        if os.path.isdir(src):
            print(f"[INFO] Moving folder {src} → {dst}")
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.move(src, dst)
        else:
            print(f"[WARN] Folder {src} not found, skipping.")

        # Step 2: For this targetfolder, check for a PDB inside outfiles/{targetfolder}
        if os.path.isdir(dst):
            pdb_found = False
            for file in os.listdir(dst):
                if file.lower().endswith(".pdb"):
                    src_pdb = os.path.join(dst, file)
                    dst_pdb = os.path.join(mhc_pdbs_dir, f"{targetfolder}.pdb")
                    print(f"[INFO] Moving PDB {src_pdb} → {dst_pdb}")
                    shutil.move(src_pdb, dst_pdb)
                    pdb_found = True
                    break  # only use the first .pdb per targetfolder
            if not pdb_found:
                print(f"[WARN] No PDB found in {dst}")

if __name__ == "__main__":
    main()

