function [mmC, mC, rpfC, rC, rbC] = getData(datadir)

    % tracks number of each molecule (indexed) at each TIME FRAME
    f1 = load([datadir, '/newmolhistperframe_1.dat']);
    mC = spconvert(f1);
%     mC = full(spconvert(f1));
    
    % tracks number of each reaction (indexed) at each TIME FRAME
    f2 = load([datadir, '/newreactperframe_1.dat']);
    rpfC = spconvert(f2);
%     rpfC = full(spconvert(f2));
     
    % tracks stoichiometry of each molecule (indexed) for each REACTION
    f3 = load([datadir, '/newreactbasis_all.dat']);
    rC = spconvert(f3);
%     rC = full(spconvert(f3));
    
    % tracks stoichiometry of each reactant (indexed) for each REACTION
    f4 = load([datadir, '/newreactantsbasis_all.dat']);
    rbC = spconvert(f4);
%     rbC = full(spconvert(f4));

%     bD = load([datadir, '/bondsdict_all.txt']);
    f5 = load([datadir, '/molmolhistperframe_1.dat']);
    mmC = spconvert(f5);
    

end