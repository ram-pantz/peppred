##AFFT-HLA3DB + Similarity Predictor ##


## SETUP ##
This engine requires the use of multiple conda environments as well as a set up script to be fully operational

1. Create the 3 environments from the yml files in the main directory. Other required dependencies / prerequisites include Rorsetta and "Tensorrt"

2. After creating the conda environments fill in the environment variables in setup.py and run

3. Download params file from https://onedrive.live.com/?id=ED57FA73D19DB823%2112884&resid=ED57FA73D19DB823%2112884&e=BPt1QX&migratedtospo=true&redeem=aHR0cHM6Ly8xZHJ2Lm1zL3UvcyFBaU80bmRGei1sZnQ1RlR6MVNLRmdmSS1mNEN4P2U9QlB0MVFY&cid=ed57fa73d19db823&v=validatepermission and populate it in the params folder in AFFT-HLA3DB

## USAGE ##

#STEP 1: STRUCTURE PREDICTION IN THE AFFT-HLA3DB FOLDER ##
1. Populate the `input_seq` folder with the HLA and peptide sequence (on separate lines, see example in `input_seq/7P4B_seq.txt`). Make sure the file name is `XYZ_seq.txt` where XYZ is the target name.
2. `bash predict_structure.sh` - will run based off of files in `input_seq`. This should take 3-5 minutes on a GPU if requesting 64 GB of memory. The script is currently set up to run on a slurm scheduler and may need to be modified for your setup.
3. Outputs will be saved as a pdb to a folder named for the input. These folders will be automatically moved to outfile when dscore.sh is run

#STEP 2: PAIR GENERATION IN THE COMPARE FOLDER #
1. Run Compare/dscore.sh

