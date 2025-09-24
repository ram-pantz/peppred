#!/bin/bash
#SBATCH --job-name=genlist
#SBATCH --array=0-4
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=04:00:00
#SBATCH --output=logs/genlist_%A_%a.out
#SBATCH --error=logs/genlist_%A_%a.err

mkdir -p outputs/genlist tmp_fastas
CHUNK=chunks/pdb_chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID)

ENV=$(python -c "from misc import constants; print(constants.condapath)")
echo "From Python: $ENV"
source $ENV
conda activate compare


echo "========== GENLIST START (Task $SLURM_ARRAY_TASK_ID) =========="

# 1. Convert PDBs in this chunk to fasta
TMP_FASTA=tmp_fastas/genlist_chunk_${SLURM_ARRAY_TASK_ID}.fasta
python genfasta.py "$CHUNK" "$TMP_FASTA"

# 2. Run genlist on that fasta to produce partial CSV
OUT_CSV=outputs/genlist/genlist_part$(printf "%02d" $SLURM_ARRAY_TASK_ID).csv
python genlist.py "$TMP_FASTA" "$OUT_CSV"

echo "[TWO] Task $SLURM_ARRAY_TASK_ID wrote → $OUT_CSV"

# 3. Merge all partial CSVs into MHCs.csv (only once, after all tasks finish)
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    echo "[TWO] Waiting for other tasks to finish..."
    scontrol show job $SLURM_JOB_ID
    sleep 30  # buffer to let others complete

    echo "[TWO] Merging CSVs..."
    head -n 1 outputs/genlist/genlist_part00.csv > MHC_pdbs/MHCs.csv
    tail -n +2 -q outputs/genlist/genlist_part*.csv >> MHC_pdbs/MHCs.csv
    echo "✅ Created MHC_pdbs/MHCs.csv"
fi

echo "========== GENLIST END (Task $SLURM_ARRAY_TASK_ID) =========="
