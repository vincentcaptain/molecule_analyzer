function result = is_in_cycle(atom_1, atom_2, bond, config)
result = 0;
neighbors_1 = find((bond(1,atom_1,:)>0));
% bond(1,atom_1,:)
% atom_2
res = zeros(size(neighbors_1));
for i = 1:size(res,1)
   if bond(1,atom_1,neighbors_1(i))==atom_2
       
       res(i) = 1;
   else
       dist = getAtomsDist(bond(1,atom_1,neighbors_1(i)), atom_2, atom_1, 0, bond);
       if dist == 0
            res(i) = 0;
       else
           res(i)=dist +1;
       end
   end
end
res
[countH, countC,countCC, countCH, countHH,  visited] = createMolecule(atom_1, 0, 1,0, 0, 0, 0, 0,zeros(size(bond,2)), config, bond);
mol = [countH, countC,countCC, countCH, countHH]
        
if size(find(res>0),1)>1
    result = 1;
end
end