function [mol_name_new, molperframe] = getMoleculeFirstFrame(bondfinal, configfinal, mol_name)
% Get the list of all the molecules encountered in mol_name_new and their
% number in molperframe for the first frame.

num_atoms = size(bondfinal,1);



mol_name_new = mol_name;
mol_count = size(mol_name_new,1);
molperframe = zeros(1,mol_count);


count_molecule = 1;
molecules = strings(1,1);
visited = zeros(num_atoms);
for j = 1:num_atoms
    countH = 0;
    countC = 0;
    countCC = 0;
    countCH = 0;
    countHH = 0;
    [countH, countC,countCC, countCH, countHH,  visited] = createMolecule(j, 0, countH, countC, countCC, countCH, countHH, visited, configfinal, bondfinal);
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
        molecules(count_molecule) = st;
        count_molecule = count_molecule +1;
    end
end
[C,ia,ic] = unique(molecules);
a_counts = accumarray(ic,1);

for j = 1:size(C,2)
if ~ismember(mol_name_new, C(1,j))
    mol_count = mol_count +1;
    mol_name_new(mol_count,1) = C(1,j);
    molperframe(1, mol_count) = zeros(1,1);
end
for k =1:size(mol_name_new,1)
    if mol_name_new(k,1) == C(1,j)
        mol_number = k;
    end
end
molperframe(1,mol_number) = a_counts(j);
end
end

