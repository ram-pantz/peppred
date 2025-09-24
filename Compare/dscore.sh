#!/bin/bash
#SBATCH --job-name=mhc_driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/driver.out
#SBATCH --error=logs/driver.err

mkdir -p logs outputs chunks
mkdir -p storedfiles
rm logs/*
rm outputs/genfasta/*
rm outputs/genlist/*
rm chunks/*
rm tmp_fastas/*
rm peptide_dihedrals/*.csv
rm anchor_class.csv

echo "========== DRIVER START =========="
# CONDA ACTIVATIONS #
ENV=$(python -c "from misc import constants; print(constants.condapath)")
echo "From Python: $ENV"
source $ENV
conda activate compare

# SET UP DRIVER #
LOC=$(python -c "from misc import constants; print(constants.mainpath)")
echo "From Python: $LOC"
python $LOC/AFFT-HLA3DB/store.py $LOC
python $LOC/AFFT-HLA3DB/npz.py $LOC

if [ -d "$LOC/Compare/MHC_pdbs" ]; then
    mv "$LOC/Compare/MHC_pdbs" "$LOC/Compare/storedfiles"
else
    echo "[WARN] no previous MHCs dir found in Compare directory"
fi

if [ -d "$LOC/Compare/scorings" ]; then
    mv "$LOC/Compare/scorings" "$LOC/Compare/storedfiles"
else
    echo "[WARN] no previous NPZ dir found in Compare directory"
fi

if [ -d "$LOC/AFFT-HLA3DB/MHC_pdbs" ]; then
    mv "$LOC/AFFT-HLA3DB/MHC_pdbs" "$LOC/Compare/"
else
    echo "[WARN] Failed to find MHC_pds Dir in AFFT-HLA3DB"
fi

if [ -d "$LOC/AFFT-HLA3DB/scorings" ]; then
    mv "$LOC/AFFT-HLA3DB/scorings" "$LOC/Compare/"
else
    echo "[WARN] Failed to find Scorings Dir in AFFT-HLA3DB"
fi

## BEGIN DRIVER ##
echo "Reordering PDBs"
python reorder.py
echo "Reordering Done"

echo "[DRIVER] Splitting PDB list..."
ls MHC_pdbs/*.pdb > chunks/pdb_list.txt
split -n l/5 --numeric-suffixes=0 --suffix-length=2 chunks/pdb_list.txt chunks/pdb_chunk_
echo "[DRIVER] Created pdb chunks"

# STEP 1: Submit genfasta array
jid1=$(sbatch --parsable --array=0-4 one.sh)
echo "[DRIVER] Submitted genfasta job array: $jid1"

# STEP 2: Submit genlist array, dependent on genfasta
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --array=0-4 two.sh)
echo "[DRIVER] Submitted genlist job array: $jid2"

# STEP 3: Submit rest pipeline, dependent on genlist
jid3=$(sbatch --parsable --dependency=afterok:$jid2 three.sh)
echo "[DRIVER] Submitted rest pipeline job: $jid3"

echo "========== DRIVER END =========="

