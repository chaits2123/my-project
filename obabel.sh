for f in ligand*.pdb; do
    obabel "$f" -O "${f%.pdb}.pdbqt" --partialcharge gasteiger
done

