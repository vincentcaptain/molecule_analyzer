function result = isPossible(mu, config, bond)
result = 0;
global reactset;
if sum(reactset(mu,:)) >0
    result = 1;
    return;
else
    [x, store] = getNumOfMol(bond, config);
    mol = find(reactset(mu,:)<0);
    if reactset(mu,mol(1)) == -2
        for i = 1:size(find(store(mol(1),:,1)>0),2)
            for j = i+1:size(find(store(mol(1),:,1)>0),2)
                if ~ismember(store(mol(1), j,1),bond(1,store(mol(1),i,1),:)) 
                    result = 1;
                    return
                end
            end
        end
    else
        for i = 1:size(find(store(mol(1),:,1)>0),2)
            for j = 1:size(find(store(mol(2),:,1)>0),2)
                if ~ismember(store(mol(2), j,1),bond(1,store(mol(1),i,1),:)) 
                    result = 1;
                    return
                end
            end
        end
    end
end


end