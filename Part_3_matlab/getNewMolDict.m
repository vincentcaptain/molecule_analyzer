function newmoldict = getNewMolDict(datadir)


newmoldict = zeros(1,2);
fileID = fopen([datadir, '/newmoldict.txt']);
moleculeline = fgets(fileID);
molid = 1;
while ischar(moleculeline)  %% stuff exists
        if size(strfind(moleculeline, '('))>0    %% if molname appears in the line
            l = strsplit(moleculeline, ')');
            l = strsplit(string(l(1)), '(');
            mol1 = strsplit(string(l(2)), ' ');
            mol1 = str2double(mol1(1));
            mol2 = strsplit(string(l(3)), ' ');
            mol2 = str2double(mol2(1));
            newmoldict(molid, 1) = mol1;
            newmoldict(molid,2)= mol2;
        else
            l = strsplit(moleculeline, ' ');
            mol = str2double(l(3));
            newmoldict(molid, 1) = mol;
        end
        moleculeline = fgets(fileID);
        molid = molid + 1;
end
fclose(fileID);

end