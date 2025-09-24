#!/bin/bash
#SBATCH --job-name=genfasta
#SBATCH --array=0-4
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=04:00:00
#SBATCH --output=logs/genfasta_%A_%a.out
#SBATCH --error=logs/genfasta_%A_%a.err

ENV=$(python -c "from misc import constants; print(constants.condapath)")
echo "From Python: $ENV"
source $ENV
conda activate compare

mkdir -p outputs/genfasta tmp_fastas

CHUNK=chunks/pdb_chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID)

echo "[ONE] Task $SLURM_ARRAY_TASK_ID processing $CHUNK"
python genfasta.py "$CHUNK" "outputs/genfasta/genfasta_part$(printf "%02d" $SLURM_ARRAY_TASK_ID).fasta"


if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    echo "[ONE] Waiting for other tasks..."
    scontrol show job $SLURM_JOB_ID
    sleep 30
    cat outputs/genfasta/genfasta_part*.fasta > MHC_pdbs/MHCs.fasta
    echo "✅ Created MHC_pdbs/MHCs.fasta"
fi

