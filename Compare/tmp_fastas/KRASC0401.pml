
reinitialize
load MHC_pdbs/KRASC0401_reordered.pdb, mol
if cmd.count_atoms("mol and chain A") > 0:
    save tmp_fastas/KRASC0401_A.fasta, mol and chain A
if cmd.count_atoms("mol and chain B") > 0:
    save tmp_fastas/KRASC0401_B.fasta, mol and chain B
quit
