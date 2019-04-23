function result = is_in_cycle(atom_1, atom_2, bond, config)
% Check if atom_1 and atom_2 are in a cycle (you know that they are in the
% same molecule).

result = 0;
neighbors_1 = find((bond(atom_1,:)>0));
res = zeros(size(neighbors_1));
for i = 1:size(res,1)
   if bond(atom_1,neighbors_1(i))==atom_2
       
       res(i) = 1;
   else
       dist = getAtomsDist(bond(atom_1,neighbors_1(i)), atom_2, atom_1, 0, bond);
       if dist == 0
            res(i) = 0;
       else
           res(i)=dist +1;
       end
   end
end
[countH, countC,countCC, countCH, countHH,  visited] = createMolecule(atom_1, 0, 1,0, 0, 0, 0, 0,zeros(size(bond,1)), config, bond);
mol = [countH, countC,countCC, countCH, countHH];
        
if size(find(res>0),1)>1
    result = 1;
end
end