function atoms_in_mol = getAtomInMol(atom,bond, atoms_in_mol)
if ~ismember(atom,atoms_in_mol)
    atoms_in_mol(size(atoms_in_mol,1)+1,1) = atom;
    for i = 1:size(find(bond(1,atom,:)>0))
        atoms_in_mol = getAtomInMol(bond(1,atom,i), bond, atoms_in_mol);
    end
end