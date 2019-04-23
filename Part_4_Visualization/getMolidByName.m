function molid = getMolidByName(molname, datadir, tag)

    % PULL FROM TEXT FILE
    if tag == 0
        fileID = fopen([datadir, '/moleculedict_all.txt']);
    end
    if tag == 1
        fileID = fopen([datadir, '/molmoleculedict_all.txt']);
    end
    moleculeline = fgets(fileID);   %% reads line by line
    molid = 1;
    while ischar(moleculeline)  %% stuff exists
        if strfind(moleculeline, molname)    %% if molname appears in the line
            break;
        end
        moleculeline = fgets(fileID);
        molid = molid + 1;
    end
end