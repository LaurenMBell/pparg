#the bash script for running AutoDock Vina!! 

# preparing the receptor
# making a PDBQT file for the target protein
mk_prepare_receptor.py -i 1iep_receptorH.pdb -o 1iep_receptor -p -v \ 
    --box_size 20 20 20 --box_center 15.190 53.903 16.917

# preparing the ligand
# making a PBDQT file for the input ligand molecule
mk_prepare_ligand.py -i 1iep_ligand.sdf -o 1iep_ligand.pdbqt

# running autodock vina 
# NEED: 