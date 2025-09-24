#       Sgourakis Lab
#   Author: Ram Pantula
#   Date: July 1, 2025
#   Email: rpantula@sas.upenn.edu

#!/usr/bin/env python3
import sys, csv, re, os

if len(sys.argv) != 3:
    print("Usage: python genlist.py <fasta_chunk> <output_csv>")
    sys.exit(1)

fasta_path = sys.argv[1]
out_csv = sys.argv[2]
imgt_fasta_path = "imgt_allelic_freq_0.05_NMDP.fasta"

os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)

def read_fasta(path):
    records = []
    if not os.path.exists(path):
        return records
    with open(path, "r") as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records

# Build allele DB: exact sequence -> allele name from IMGT file
allele_db = {}
for header, seq in read_fasta(imgt_fasta_path):
    allele_db[seq] = header

# Parse chunk fasta and group chains by PDB
group = {}  # pdb -> {"A": seq, "B": seq}
for header, seq in read_fasta(fasta_path):
    pdb = None; chain = None
    # Preferred: "<PDB> Chain X"
    m = re.match(r"^(.+?)\s+Chain\s+([AB])$", header)
    if m:
        pdb, chain = m.group(1), m.group(2)
    else:
        # Fallback: "<PDB>_X"
        m = re.match(r"^(.+?)_([AB])$", header)
        if m:
            pdb, chain = m.group(1), m.group(2)
        else:
            # Last resort: strip a trailing _A/_B or space-A/B
            if header and header[-1] in "AB":
                pdb = re.sub(r"[_\s][AB]$", "", header)
                chain = header[-1]
    if not pdb or not chain:
        continue
    group.setdefault(pdb, {})[chain] = seq

rows = []
for pdb, chains in group.items():
    if "A" not in chains or "B" not in chains:
        continue
    allele = allele_db.get(chains["A"], "UNKNOWN")
    pep = chains["B"]
    rows.append([pdb, allele, len(pep), pep])

with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["PDB", "Allele", "Peptide length", "Peptide Sequence"])
    w.writerows(rows)

print(f"Wrote {len(rows)} rows → {out_csv}")

