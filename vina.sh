#!/bin/bash

# AutoDock Vina Batch Docking Script

PROTEIN="protein.pdbqt"

CENTER_X=10.0
CENTER_Y=25.0
CENTER_Z=35.0
SIZE_X=20.0
SIZE_Y=20.0
SIZE_Z=20.0

echo "Converting ligands to PDBQT..."
for f in ligand*.pdb; do
    echo "Converting $f..."
    obabel "$f" -O "${f%.pdb}.pdbqt" --partialcharge gasteiger
done

echo "Starting Vina docking..."
for LIGAND in ligand*.pdbqt; do
    BASENAME=$(basename $LIGAND .pdbqt)
    echo "Docking $LIGAND ..."
    vina --receptor $PROTEIN \
         --ligand $LIGAND \
         --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
         --size_x $SIZE_X --size_y $SIZE_Y --size_z $SIZE_Z \
         --out ${BASENAME}_out.pdbqt \
         --log ${BASENAME}.log
done

echo "Docking complete!"

