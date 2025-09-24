#       Sgourakis Lab
#   Author: Ram Pantula
#   Date: July 1, 2025
#   Email: rpantula@sas.upenn.edu

#!/usr/bin/env python3
import os
import sys
import numpy as np

# --------------------------------------------------------------------
# CONFIG
if len(sys.argv) > 1:
    BASE_DIR = os.path.abspath(sys.argv[1])
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_ROOT = os.path.join(BASE_DIR, "AFFT-HLA3DB/outfiles")
STORAGE_ROOT = os.path.join(BASE_DIR, "AFFT-HLA3DB/scorings")

print(f"[INFO] Base directory set to: {BASE_DIR}")
print(f"[INFO] Input folder: {INPUT_ROOT}")
print(f"[INFO] Output folder: {STORAGE_ROOT}")
# --------------------------------------------------------------------

def combine_npy_to_npz(dirpath, npy_files, out_path):
    arrays = {}
    for fname in npy_files:
        key = os.path.splitext(fname)[0]
        fpath = os.path.join(dirpath, fname)
        try:
            arrays[key] = np.load(fpath)
        except Exception as e:
            print(f"[!] Failed to load {fpath!r}: {e}")

    if arrays:
        np.savez_compressed(out_path, **arrays)
        print(f"[+] Wrote {out_path}")
    else:
        print(f"[!] No valid arrays to save in {dirpath!r}")

def main():
    input_root = os.path.abspath(INPUT_ROOT)
    storage_root = os.path.abspath(STORAGE_ROOT)
    os.makedirs(storage_root, exist_ok=True)

    print(f"Scanning for .npy files under {input_root!r} …")
    for dirpath, _, files in os.walk(input_root):
        # find any .npy files in this folder
        npy_files = [f for f in files if f.lower().endswith(".npy")]
        if not npy_files:
            continue

        folder_name = os.path.basename(dirpath)
        out_filename = f"{folder_name}.npz"
        out_path = os.path.join(storage_root, out_filename)

        print(f"Combining {len(npy_files)} .npy files in {dirpath!r} → {out_path!r}")
        combine_npy_to_npz(dirpath, npy_files, out_path)

if __name__ == "__main__":
    main()

