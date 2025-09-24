##AFFT-HLA3DB + Similarity Predictor ##


## SETUP ##
This engine requires the use of multiple conda environments as well as a set up script to be fully operational
Create the 3 environments from the yml files below
Other required dependencies / prerequisites include Rorsetta and "Tensorrt"

After creating the conda environments fill in the environment variables in setup.py and run

## USAGE ##

#STEP 1: STRUCTURE PREDICTION IN THE AFFT-HLA3DB FOLDER ##
1. Populate the `input_seq` folder with the HLA and peptide sequence (on separate lines, see example in `input_seq/7P4B_seq.txt`). Make sure the file name is `XYZ_seq.txt` where XYZ is the target name.
2. `bash predict_structure.sh` - will run based off of files in `input_seq`. This should take 3-5 minutes on a GPU if requesting 64 GB of memory. The script is currently set up to run on a slurm scheduler and may need to be modified for your setup.
3. The output can be found in a folder named `XYZ` where XYZ is the target name. See example in `7P4B/` where `outfile_7P4B_model_1_model_2_ptm_ft.pdb` is the final structure model.

#STEP 2: PAIR GENERATION IN THE COMPARE FOLDER #
1. Ensure all conda environments are installed
2. Run Compare/dscore.sh

