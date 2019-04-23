function [x_new,x0, bond_new, config_new, mol_name_new, store_new] = makeReaction(bond, config, mu, mol_name, x_mol, x0_old, store_old)
% Make the reaction happen and update all the variables.

store = store_old;
x = x0_old;
global reactset;
config_new = config;
bond_new = bond;
if sum(reactset(mu,:)) <0
    mol = find(reactset(mu,:)<0);
    atom = zeros(2,1);
    for i = 1:size(mol,2)
        if reactset(mu,mol(i)) == -2
            atom(1,1) = store(mol(i), ceil(x(mol(i),1)*rand),1);
            atom(2,1) = store(mol(i), ceil(x(mol(i),1)*rand),1);
            while atom(1,1) == atom(2,1)
                atom(2,1) = store(mol(i), ceil(x(mol(i),1)*rand),1);
            end
            test = 0;
            counter_1 = 0;
            while test==0
                j = 1;
                counter_2 = 0;
                while test == 0 && j <= size(find(store(mol(1), :,1)>0),2)
                    k = ismember(store(mol(1), j,1),bond(atom(1),:));
                    if ~k
                        test = 1;
                    end
                    j = j+1;
                    counter_2 = counter_2 +1;
                    if counter_2>1000
                        disp('BIG ERROR')
                        break
                    end
                end
                if test == 0
                    atom(1) = store(mol(1), ceil(x(mol(1),1)*rand),1);
                end
                if counter_1>1000
                    disp('OTHER BIG ERROR')
                    break
                end
            end
            while ismember(atom(2), bond(atom(1),:))
                atom(2) = store(mol(1), ceil(x(mol(1),1)*rand),1);
            end    
        else
            atom(i) = store(mol(i), ceil(x(mol(i),1)*rand),1);
            if i==1
                test = 0;
                counter_1 = 0;
                while test==0
                    j = 1;
                    counter_2 = 0;
                    while test == 0 && j <= size(find(store(mol(2), :,1)>0),2)
                        k = ismember(store(mol(2), j,1),bond(atom(1),:));
                        if ~k
                            test = 1;
                        end
                        j = j+1;
                        counter_2 = counter_2 +1;
                        if counter_2>1000
                            disp('BIG ERROR')
                            break
                        end
                    end
                    if test == 0
                        atom(i) = store(mol(i), ceil(x(mol(i),1)*rand),1);
                    end
                    if counter_1>1000
                        disp('OTHER BIG ERROR')
                        break
                    end
                end
            end
            if i ==2
                while ismember(atom(2), bond(atom(1),:))
                    atom(i) = store(mol(i), ceil(x(mol(i),1)*rand),1);
                end
                
            end
        end
    end
    if ismember(atom(2), bond(atom(1),:))
        disp('Careful');
    end
    if config_new(atom(1))>= 2000 && config_new(atom(2))>= 2000
        modif = 10;
    elseif config_new(atom(1))< 2000 && config_new(atom(2))<2000
        modif = 100;
    else
        modif = 1;
    end
    for i = 1:2
        j = 1;
        while bond_new(atom(i),j)>0
            j = j+1;
        end
        bond_new(atom(i),j) = atom(3-i);
        config_new(atom(i)) = config_new(atom(i)) + modif;
    end
else
    mol = find(reactset(mu,:)<0);
    atom = zeros(2,1);
    k = ceil(x(mol,1)*rand);
    atom(1) = store(mol, k,1);
    atom(2) = store(mol, k,2);
    if config_new(atom(1))>= 2000 && config_new(atom(2))>= 2000
        modif = 10;
    elseif config_new(atom(1))< 2000 && config_new(atom(2))<2000
        modif = 100;
    else
        modif = 1;
    end
    for i = 1:2
        j = 1;
        reax = 0;
        while bond_new(atom(i),j)>0
            if reax == 1
                bond_new(atom(i),j-1) = bond_new(atom(i),j);
                bond_new(atom(i),j) = 0;
            end
            if bond_new(atom(i),j) == atom(3-i)
                bond_new(atom(i),j) = 0;
                reax = 1;
            end
            j = j+1;
         end
        config_new(atom(i)) = config_new(atom(i)) - modif;
    end
end

[mol_name_new, x_new] = getMoleculePerFrame(bond_new, config_new, mol_name, bond, config, x_mol, atom);
[x0_2, store_new_2] = getFirstNumOfMol(bond_new, config_new);
[x0, store_new] = getNumOfMol(bond_new, config_new, bond, config, x, store, atom);
for i=1:size(x0,1)
    if x0(i) ~= x0_2(i)
        disp(i)
    end
end


end


