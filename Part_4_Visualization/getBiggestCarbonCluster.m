function [biggest_CC, num_CC, numofCinCC] = getBiggestCarbonCluster(mD, datadir, tag, molecule_name)
    if tag == 1
        name = getNameByMolid(1:size(mD,2),datadir);
        num_of_mol = size(name,2);
    else
        name = molecule_name;
        num_of_mol = size(name,1);
    end
    num_C = zeros(num_of_mol,1);
    biggest_CC = zeros(size(mD,1),1);
    num_CC = zeros(size(mD,1),1);
    numofCinCC = zeros(size(mD,1),1);
    for i = 1:num_of_mol
        n = name(i);
        n = n{:};
        n = strsplit(n, ' ');
        n = n(1);
        n = n{:};
        if size(n,1)>0
            if n(1) == 'C'
                num_C(i) = str2num(erase(n,'C'));
            end
        end
    end    

    for i = 1:size(mD,1)
        p = mD(i,:);
        pp = find(p);
        for j = 1:size(pp,2)
%             biggest_CC(i)
%             pp(j)
%             size(num_C)
%             num_C(pp(j))
            biggest_CC(i) = max(biggest_CC(i),num_C(pp(j)));
            if num_C(pp(j)) >4
                num_CC(i) = num_CC(i) +1;
                numofCinCC(i) = numofCinCC(i) + num_C(pp(j));
            end         
        end
    end