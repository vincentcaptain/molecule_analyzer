function [x, store] = getNumOfMol(bond,config)

global newmoldict;

nummol = size(newmoldict,1);

x = zeros(nummol,1);
store = zeros(1,1,2);
% fileID = fopen([datadir, '/newmoldict.txt']);
% moleculeline = fgets(fileID);
molid = 1;

for k = 1:nummol
    count = 0;
    store(molid, :, :) = 0;
    if newmoldict(k,2) ~=0
        mol1 = newmoldict(k,1);
        mol2 = newmoldict(k,2);
        for i = 1:size(config,2)
            if config(1,i) == mol1
                j = 1;
                while bond(1,i,j)>0
                    if config(1, bond(1,i,j)) == mol2
                        count = count + 1;
                        store (molid,count, 1) = i;
                        store(molid, count, 2) = bond(1, i,j);
                    end
                    j = j+1;
                end

            end
         end
    else
        mol = newmoldict(k,1);
        for i =1:size(config,2)
            if config(1,i) == mol
                count = count + 1;
                store(molid, count, 1) = i;
            end
        end
    end
    x(molid,1) = count;
    molid = molid + 1;
end


% while ischar(moleculeline)  %% stuff exists
%         count = 0;
%         store(molid, :, :) = 0;
%         if size(strfind(moleculeline, '('))>0    %% if molname appears in the line
%             l = strsplit(moleculeline, ')');
%             l = strsplit(string(l(1)), '(');
%             mol1 = strsplit(string(l(2)), ' ');
%             mol1 = str2double(mol1(1));
%             mol2 = strsplit(string(l(3)), ' ');
%             mol2 = str2double(mol2(1));
%             for i = 1:size(config,2)
%                 if config(1,i) == mol1
%                     j = 1;
%                     while bond(1,i,j)>0
%                         if config(1, bond(1,i,j)) == mol2
%                             count = count + 1;
%                             store (molid,count, 1) = i;
%                             store(molid, count, 2) = bond(1, i,j);
%                         end
%                         j = j+1;
%                     end
%                     
%                 end
%             end
%         else
%             l = strsplit(moleculeline, ' ');
%             mol = str2double(l(3));
%             for i =1:size(config,2)
%                 if config(1,i) == mol
%                     count = count + 1;
%                     store(molid, count, 1) = i;
%                 end
%             end
%         end
%         x(molid,1) = count;
%         moleculeline = fgets(fileID);
%         molid = molid + 1;
% end
% fclose(fileID);
end

