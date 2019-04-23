function [x0, store_new] = getNumOfMol(bond_new, config_new, bond, config, x0_old, store_old, atom)
%Mofidy x0 and store and update them after the reaction happened.

x0 = x0_old;
store_new = store_old;

global newmoldict;
nummol = size(newmoldict,1);

for k =1:nummol
    if newmoldict(k,2) ~=0
        mol1 = newmoldict(k,1);
        mol2 = newmoldict(k,2);
        if config(atom(1)) == mol1 && config(atom(2)) == mol2 && ismember(atom(1), bond(atom(2),:))
            x0(k,1) = x0(k,1) -1;
            store_new = removeFromStore(store_new, k, atom(1), atom(2), config);
        elseif config(atom(1)) == mol2 && config(atom(2)) == mol1 && ismember(atom(1), bond(atom(2),:))
            x0(k,1) = x0(k,1) -1;
            store_new = removeFromStore(store_new, k, atom(2), atom(1), config);
        end
        for i =1:size(atom,1)
            if config(atom(i)) == mol1
                for j=1:size(bond,2)
                    if bond(atom(i),j)>0
                        if config(bond(atom(i),j)) == mol2
                            if ~(bond(atom(i),j) == atom(1) || bond(atom(i),j) == atom(2))
                                x0(k,1) = x0(k,1) -1;
                                store_new = removeFromStore(store_new, k, atom(i), bond(atom(i),j), config);
                            end
                        end
                    end
                end
            elseif config(atom(i)) == mol2
                for j=1:size(bond,2)
                    if bond(atom(i),j)>0
                        if config(bond(atom(i),j)) == mol1
                            if ~(bond(atom(i),j) == atom(1) || bond(atom(i),j) == atom(2))
                                x0(k,1) = x0(k,1) -1;
                                store_new = removeFromStore(store_new, k, bond(atom(i),j),atom(i), config);
                            end
                        end
                    end
                end
            end
        end
    else 
        mol1 = newmoldict(k,1);
        for i =1:size(atom,1)
            if config(atom(i)) == mol1
                x0(k,1) = x0(k,1) -1;
                store_new = removeFromStore(store_new, k, atom(i), 0, config);
            end
        end
    end
end

% disp('Add')

for k =1:nummol
    if newmoldict(k,2) ~=0
        mol1 = newmoldict(k,1);
        mol2 = newmoldict(k,2);
        if config_new(atom(1)) == mol1 && config_new(atom(2)) == mol2 && ismember(atom(1), bond_new(atom(2),:))
            x0(k,1) = x0(k,1) +1;
            store_new = addToStore(store_new, k, atom(1),atom(2));
        elseif config_new(atom(1)) == mol2 && config_new(atom(2)) == mol1 && ismember(atom(1), bond_new(atom(2),:))
            x0(k,1) = x0(k,1) +1;
            store_new = addToStore(store_new, k, atom(2),atom(1));
        end
        for i =1:size(atom,1)
            if config_new(atom(i)) == mol1
                for j=1:size(bond_new,2)
                    if bond_new(atom(i),j)>0
                        if config_new(bond_new(atom(i),j)) == mol2
                            if ~(bond_new(atom(i),j) == atom(1) || bond_new(atom(i),j) == atom(2))
                                x0(k,1) = x0(k,1) +1;
                                store_new = addToStore(store_new, k, atom(i), bond_new(atom(i),j));
                            end
                        end
                    end
                end
            elseif config_new(atom(i)) == mol2
                for j=1:size(bond_new,2)
                    if bond_new(atom(i),j)>0
                        if config_new(bond_new(atom(i),j)) == mol1
                            if ~(bond_new(atom(i),j) == atom(1) || bond_new(atom(i),j) == atom(2))                          
                                x0(k,1) = x0(k,1) +1;
                                store_new = addToStore(store_new, k, bond_new(atom(i),j),atom(i));
                            end
                        end
                    end
                end
            end
        end
    else 
        mol1 = newmoldict(k,1);
        for i =1:size(atom,1)
            if config_new(atom(i)) == mol1
                x0(k,1) = x0(k,1) +1;
                store_new = addToStore(store_new, k, atom(i), 0);
            end
        end
    end
end
end

