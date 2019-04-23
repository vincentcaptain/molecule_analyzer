function atom_dist = getAtomsDist(atom_1,atom_2, father, atom_dist, bond)
% Outputs the distance between 2 atoms in a same molecule.
result = zeros(0,1);
for i = 1:size(find(bond(atom_1,:)>0))
    if bond(atom_1,i) == atom_2
        result(size(result,1)+1,1) = 1;
    elseif ismember(bond(atom_1,i),father)
        continue
    else
        new_father = father;
        new_father(size(new_father,1)+1,1) = atom_1;
        a = getAtomsDist(bond(atom_1,i), atom_2, new_father, atom_dist, bond);
        if a ~= 0
            result(size(result,1)+1,1) = a +1;
        end
    end
if size(result,1)>0
    atom_dist = min(result);
else
    atom_dist = 0;    
end

end