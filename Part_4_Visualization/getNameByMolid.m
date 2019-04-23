function molfullnames = getNameByMolid(molids,datadir)
    molids = reshape(molids,1,length(molids));

    % PULL FROM TEXT FILE
%     if tag ==1
     fileID = fopen([datadir,'/molmoleculedict_all.txt']);
%     else
%         fileID = fopen([datadir,'/molmoleculedict_all.txt']);
    
    molfullnames = [];
    %oldmolid = 0;
    for molid=molids
        % NOTE: TEXTSCAN READS FROM LAST POSITION IN FILE
	    %molline = textscan(fileID,'%s',1,'delimiter','\n','headerlines',molid-oldmolid-1);
        molline = textscan(fileID,'%s',1,'delimiter','\n','headerlines',molid-1);
    	mol = strsplit(char(molline{1}),' - ');
        if length(mol)<2
            molid
            molline
        end
    	molfullnames = [molfullnames,mol(2)];
        %oldmolid = molid;
        frewind(fileID);
    end
    
    fclose(fileID);
end