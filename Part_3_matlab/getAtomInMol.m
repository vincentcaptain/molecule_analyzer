function atoms_in_mol = getAtomInMol(atom,bond, atoms_in_mol)
% This function outputs the list of index of the atoms in the molecule
% containing 'atom'. By recursion.
if ~ismember(atom,atoms_in_mol)
    atoms_in_mol(size(atoms_in_mol,1)+1,1) = atom;
    for i = 1:size(find(bond(atom,:)>0))
        atoms_in_mol = getAtomInMol(bond(atom,i), bond, atoms_in_mol);
    end
end
end