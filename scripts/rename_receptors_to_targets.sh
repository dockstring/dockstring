for pdb in $( ls *_receptor.pdb ); do
    name=$( basename $pdb _receptor.pdb )
    mv "$name"_receptor.pdb "$name"_target.pdb
done

for pdbqt in $( ls *_receptor.pdbqt ); do
    name=$( basename $pdbqt _receptor.pdbqt )
    mv "$name"_receptor.pdbqt "$name"_target.pdbqt
done
