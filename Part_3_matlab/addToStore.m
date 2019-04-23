function store_new = addToStore(store_new, config_atom, atom_1, atom_2)
% This function adds new elements to the store matrix, because it needs to
% keep the non-zero values at the beginning and the zero values at the end.

idx_1 = find(store_new(config_atom, :, 1), 1, 'last');
if isempty(idx_1)
    idx_1 = 0;
elseif size(idx_1,2)>1
    disp('Bad storage in adding')
end

store_new(config_atom,idx_1 +1 ,1) = atom_1;
store_new(config_atom,idx_1 + 1,2) = atom_2;
end