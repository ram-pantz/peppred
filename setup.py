#!/usr/bin/env python3
import os
import re

# =============================
# User Configuration (edit here)
# =============================

# Paths
ROOT = "ENTER PATH TO PairFold DIRECTORY"
ROSETTA = "ENTER PATH TO ROSETTA MAIN FOLDER"
CONDA = "ENTER PATH TO CONDA BIN/ACTIVATE FOLDER"
TENSOR = "ENTER PATH TO TENSORRT PACKAGE"

# List of target scripts to update
TARGET_FILES = [
    f"{ROOT}/AFFT-HLA3DB/fold.sh",
    "{ROOT}/AFFT-HLA3DB/predict_structure.sh",
    "{ROOT}/Compare/misc/constant.py
]

# =============================
# Implementation
# =============================

# Collect config vars automatically (all UPPERCASE)
config_vars = {k: str(v) for k, v in globals().items() if k.isupper() and k not in ["TARGET_FILES"]}

def replace_placeholders(content):
    """Replace {{VAR}} placeholders with values from config_vars."""
    def repl(match):
        var_name = match.group(1)
        if var_name in config_vars:
            return config_vars[var_name]
        else:
            print(f"[!] Warning: {var_name} not defined in config section")
            return match.group(0)  # leave placeholder unchanged
    return re.sub(r"\{\{(\w+)\}\}", repl, content)

def main():
    for path in TARGET_FILES:
        if not os.path.exists(path):
            print(f"[!] Skipping {path} (not found)")
            continue

        with open(path, "r") as f:
            original = f.read()

        updated = replace_placeholders(original)

        # Backup before overwrite
        backup_path = path + ".bak"
        if not os.path.exists(backup_path):
            with open(backup_path, "w") as f:
                f.write(original)

        with open(path, "w") as f:
            f.write(updated)

        print(f"[+] Updated {path}")

if __name__ == "__main__":
    main()

