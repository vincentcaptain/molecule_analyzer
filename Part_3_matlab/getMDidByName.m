function molid = getMDidByName(molname, datadir)

    % PULL FROM TEXT FILE
    fileID = fopen([datadir, '/moleculedict_all.txt']);
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