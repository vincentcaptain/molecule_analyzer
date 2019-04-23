function [num_cycles, numofCincycles, numofHincycles, diff] = getCycles(mD, datadir, tag, molecule_name)
    if tag == 1
        name = getNameByMolid(1:size(mD,2),datadir);
        num_of_mol = size(name,2);
    else
        name = molecule_name;
        num_of_mol = size(name,1);
    end
    num_C = zeros(num_of_mol,1);
    num_CCbond = zeros(num_of_mol,1);
    num_H =  zeros(num_of_mol,1);
    num_cycles = zeros(size(mD,1),1);
    numofCincycles = zeros(size(mD,1),1);
    numofHincycles = zeros(size(mD,1),1);
    size(num_cycles);
    for i = 1:num_of_mol
        if size(name(i),1)>0
            n = name(i);
            n = n{:};
            n = strsplit(n, ' ');
            n = n(1, 1:end-1);
            for j = 1:size(n,2)
                k = n(j);
                k = k{:};
                if k(1) == 'C'
                    num_C(i) = str2double(erase(k,'C'));
                elseif k(1) == 'H'
                    num_H(i) = str2double(erase(k,'H'));
                elseif size(k,2)>5
                    if k(end-4:end) == '(C-C)'
                        num_CCbond(i) = str2double(erase(k,'(C-C)'));
                    end
                end
            end
        end
    end    
    diff = zeros(1,1);
    for i = 1:size(mD,1)
        p = mD(i,:);
        pp = find(p>0);
        for j = 1:size(pp,2)
            if num_C(pp(j)) <= num_CCbond(pp(j)) && num_C(pp(j))>0 && num_CCbond(pp(j))>0
                num_cycles(i) = num_cycles(i) + num_CCbond(pp(j)) - num_C(pp(j)) + 1;
                %num_cycles(i) = num_cycles(i) + 1;
                numofCincycles(i) = numofCincycles(i) + num_C(pp(j));
                numofHincycles(i) = numofHincycles(i) + num_H(pp(j));
            end         
        end
        diff(1) = 0;
        if i>1
            diff(i) = num_cycles(i)-num_cycles(i-1); 
        end
    end
end