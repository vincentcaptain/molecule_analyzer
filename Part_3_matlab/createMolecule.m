function [countH_new, countC_new,countCC_new, countCH_new, countHH_new, visited_new] = createMolecule(id, id_prev, time, countH, countC,countCC, countCH, countHH, visited, configfinal, bondfinal);
countH_new = countH;
countC_new = countC;
countCH_new = countCH;
countCC_new = countCC;
countHH_new = countHH;
visited_new = visited;

if visited_new(id) == 1
    if id_prev ~= 0
        if configfinal(time,id_prev) <2000 && configfinal(time,id) <2000
            countCC_new = countCC_new +1;
        elseif configfinal(time,id_prev) >=2000 && configfinal(time,id) >=2000
            countHH_new = countHH_new +1;
        else
            countCH_new = countCH_new +1;
        end
    end

elseif visited_new(id) == 0
    if configfinal(time,id) <2000
        countC_new = countC_new + 1;
    else
        countH_new = countH_new + 1;
    end
    if id_prev ~= 0
        if configfinal(time,id_prev) <2000 && configfinal(time,id) <2000
            countCC_new = countCC_new +1;
        elseif configfinal(time,id_prev) >=2000 && configfinal(time,id) >=2000
            countHH_new = countHH_new +1;
        else
            countCH_new = countCH_new +1;
        end
    end
    visited_new(id) = 1;
    for k = 1:10
        if bondfinal(time, id,k)>0 && bondfinal(time,id,k) ~= id_prev
            [countH_new, countC_new, countCC_new, countCH_new, countHH_new, visited_new] = createMolecule(bondfinal(time, id,k),id, time, countH_new, countC_new,countCC_new, countCH_new, countHH_new, visited_new, configfinal, bondfinal);
        end
    end
    visited_new(id) = -1;
end
end