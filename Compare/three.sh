#!/bin/bash
#SBATCH --job-name=mhc_rest
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --output=logs/rest.out
#SBATCH --error=logs/rest.err

echo "========== COMP PIPELINE START =========="
ENV=$(python -c "from misc import constants; print(constants.condapath)")
echo "From Python: $ENV"
source $ENV
conda activate compare

echo "[THREE] Running anchor_class..."
python -m anchor_class.anchor_class > logs/anchor.txt 2>&1 && \
echo "[THREE] anchor_class done."

echo "[THREE] Running dihedrals..."
python -m peptide_dihedrals.dihedrals > logs/dihedrals.txt 2>&1 && \
echo "[THREE] dihedrals done."

echo "[THREE] Running pairings..."
python -m pairings.pairings > logs/pairs.txt 2>&1 && \
echo "[THREE] pairings done."

conda activate train

echo "[THREE] Finding pairs"
python compare.py
echo "[THREE] pairing completed"

echo "========== COMP PIPELINE END =========="

