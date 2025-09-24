#!/bin/bash
#SBATCH --job-name=gen_inputs
#SBATCH --time=00:20:00
#SBATCH -p shortq
#SBATCH --mem=8G
#SBATCH -o inputgen.out
#SBATCH --error=inputgen.err

# To run on CHOP's Respublica Cluster, uncomment the following two lines:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{TENSOR}

source {CONDA}
conda activate alphafold

# This will likely need to be obtained from an external link
params=params/7WKJ_af_mhc_params_2351.pkl

mkdir -p slurm_logs

# Make folders and alignment files
python initialize.py

for inputfile in ./input_seq/*.txt; do
    targname=$(basename "$inputfile" | cut -f 1 -d '_')
    olog="runlogs/${targname}.out"
    elog="runlogs/${targname}.err"
    sbatch --partition=gpuq --gres=gpu:1 --mem=64G --output="$olog" --error="$elog" predict_structure.sh "$targname" "$params"
done
