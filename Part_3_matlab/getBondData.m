function [bond,config] = getBondData(datadir, numatoms, startframe);

k = importdata([datadir, '/bondperframe_1.dat']);
l = load([datadir, '/configperframe_1.dat']);
l = spconvert(l);

bond = zeros(numatoms,10);
config = zeros(numatoms,1);

for i = numatoms*(startframe-1):numatoms*startframe
    x = strsplit(string(k(i,:)),';');
    if string(x(1)) == int2str(startframe)
        a = strsplit(string(x(3)), '[');
        a = strsplit(string(a(2)), ']');
        if contains(a(1),', ')
            a = strsplit(string(a(1)), ', ');
        else
            a = a(1);
        end
        for j = 1:size(a, 2)
            bond(str2double(x(2)), j) = a(j);
        end
    end
    
    for i = 1:numatoms
        config(i) = l(startframe, i);
    
end
end