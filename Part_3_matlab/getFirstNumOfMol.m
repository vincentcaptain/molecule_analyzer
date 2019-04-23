function [x, store] = getFirstNumOfMol(bond,config)
% From bond and config, you output x which is the number of atom in each
% configuration and store which is a table where the first index is the
% index of the atom configuration of interest. Then you store the indexes
% of the atoms which are in this configuration. If the configuration is a
% bond, you store 2 atoms, if it is just one atom you store one atom and
% put 0 in the other cell.
global newmoldict;

nummol = size(newmoldict,1);

x = zeros(nummol,1);
store = zeros(1,1,2);
molid = 1;

for k = 1:nummol
    count = 0;
    store(molid, :, :) = 0;
    if newmoldict(k,2) ~=0 
        mol1 = newmoldict(k,1);
        mol2 = newmoldict(k,2);
        if mol1 ~= mol2 % when mol1 == mol2 you count twice the configuration by this technique
            for i = 1:size(config,1)
                if config(i) == mol1
                    j = 1;
                    while bond(i,j)>0
                        if config(bond(i,j)) == mol2
                            count = count + 1;
                            store(molid,count, 1) = i;
                            store(molid, count, 2) = bond(i,j);
                        end
                        j = j+1;
                    end

                end
            end
        else
            for i =1:size(config,1)-1
                for j = i+1:size(config,1)
                   if config(i) == mol1 && config(j) == mol1
                       if ismember(j, bond(i,:))
                           count = count +1;
                           store(molid,count, 1) = i;
                           store(molid, count, 2) = j;
                       end
                   end  
                end
            end
        end
    else
        mol = newmoldict(k,1);
        for i =1:size(config,1)
            if config(i) == mol
                count = count + 1;
                store(molid, count, 1) = i;
            end
        end
    end
    x(molid,1) = count;
    molid = molid + 1;
end
end

