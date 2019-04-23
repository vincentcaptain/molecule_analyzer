function [mol_name_new, molperframe] = getMoleculePerFrame(bondfinal, configfinal, mol_name, bond_old, config_old, molperframe_old, atom)
% Modify mol_name_new and molperframe to get the updated value.

num_atoms = size(bondfinal,1);
old_mols = strings(1,1);
new_mols = strings(1,1);
visited = zeros(num_atoms);
mol_name_new = mol_name;

count_molecule = 1;
for j = 1:size(atom,1)
    countH = 0;
    countC = 0;
    countCC = 0;
    countCH = 0;
    countHH = 0;
    [countH, countC,countCC, countCH, countHH,  visited] = createMolecule(atom(j), 0, countH, countC, countCC, countCH, countHH, visited, config_old, bond_old);
    count = [countC, countH,countCC, countCH, countHH];
    name = ["C", "H", "(C-C)", "(H-C)", "(H-H)"]; 
    if find(count>0)>0
        st = "";
        for k = 1:5
            if count(k)>0
                if k<=2
                    st = st + name(k) + int2str(count(k)) + " ";
                else
                    st = st + int2str(count(k))+ name(k) + " ";
                end
            end
        end
        old_mols(count_molecule) = st;
        count_molecule = count_molecule +1;
    end
end



visited = zeros(num_atoms);
count_molecule = 1;
for j = 1:size(atom,1)
    countH = 0;
    countC = 0;
    countCC = 0;
    countCH = 0;
    countHH = 0;
    [countH, countC,countCC, countCH, countHH,  visited] = createMolecule(atom(j), 0, countH, countC, countCC, countCH, countHH, visited, configfinal, bondfinal);
    count = [countC, countH,countCC, countCH, countHH];
    name = ["C", "H", "(C-C)", "(H-C)", "(H-H)"]; 
    if find(count>0)>0
        st = "";
        for k = 1:5
            if count(k)>0
                if k<=2
                    st = st + name(k) + int2str(count(k)) + " ";
                else
                    st = st + int2str(count(k))+ name(k) + " ";
                end
            end
        end
        new_mols(count_molecule) = st;
        count_molecule = count_molecule +1;
    end
end

for i=1:size(new_mols,2)
    if ~ismember(new_mols(1,i),mol_name_new)
        mol_name_new(size(mol_name_new,1)+1,1) = new_mols(1,i);
    end
end

molperframe = zeros(1,size(mol_name_new,1));
molperframe(1,1:size(molperframe_old,2)) = molperframe_old;


for i =1:size(old_mols,2)
    for j=1:size(mol_name_new,1)
        if mol_name_new(j,1) == old_mols(1,i)
            molperframe(1,j) = molperframe(1,j) -1;
        end
    end
end

for i =1:size(new_mols,2)
    for j=1:size(mol_name_new,1)
        if mol_name_new(j,1) == new_mols(1,i)
            molperframe(1,j) = molperframe(1,j) +1;
        end
    end
end
end
