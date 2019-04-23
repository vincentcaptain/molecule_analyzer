function store_new = removeFromStore(store_new, config_atom, atom_1, atom_2, config)
%Update store by removing a molecule but keeping the structure.

idx_1 = find(store_new(config_atom, :, 1) == atom_1 & store_new(config_atom,:,2) == atom_2);

if size(idx_1,2) ~= 1
    if size(idx_1,2) == 0 && config(atom_1) == config(atom_2)
        idx_1 = find(store_new(config_atom, :, 1) == atom_2 & store_new(config_atom,:,2) == atom_1);
    else
        disp('Bad store update!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
end

idx_2 = find(store_new(config_atom, :, 1), 1, 'last');

store_new(config_atom,:,:);
idx_1;
idx_2;

store_new(config_atom,idx_1,1) = store_new(config_atom,idx_2,1);
store_new(config_atom,idx_1,2) = store_new(config_atom,idx_2,2);
store_new(config_atom,idx_2,1) = 0;
store_new(config_atom,idx_2,2) = 0;

store_new(config_atom,:,:);

end