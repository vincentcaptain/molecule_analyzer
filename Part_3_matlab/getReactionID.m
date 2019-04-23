function index = getReactionID(datadir, reaction_name)
index = 0;
fileID = fopen([datadir, '/newreactdict_all.txt']);
moleculeline = fgets(fileID);   %% reads line by line
molid = 1;
while ischar(moleculeline)  %% stuff exists
    l = strsplit(moleculeline, ': ');
    l = strsplit(char(l(2)), '\r');
    if char(l(1)) == char(reaction_name)
        index = molid;
        break
    end
    
    moleculeline = fgets(fileID);
    molid = molid +1;
end


end