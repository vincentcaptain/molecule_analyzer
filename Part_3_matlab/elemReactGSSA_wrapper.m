function [tfinal, xfinal, mfinal,numt , mol_name] = elemReactGSSA_wrapper(bond0, config0,numatoms, x0, k, simtype, numframes, tscaleps)

    %% Rate constants and reaction stoichiometry
    % rate constants are input as k
    global reactset;    %% full reaction coefficients
    stoich_matrix = [
        reactset;
    ]; % rows: reactions; columns: molecules

    %% Initial Conditions
    twindow = tscaleps * min(100, numframes); % num timesteps to simulate each time, due to memory constraints
    tspan = [0, twindow]; %timesteps
    maxnumreactest = ceil(100000 * tscaleps / 0.48);  %% optimized (estimated) for memory allocation
    mol_name = strings(1,1);
    
    tfinal = zeros(ceil(numframes/10),1);
    xfinal = zeros(ceil(numframes/10),300);
    mfinal = zeros(ceil(numframes/10),1);
    [mol_name,x] = getMoleculeFirstFrame(bond0, config0, mol_name);
    [x0, store] = getFirstNumOfMol(bond0, config0);

    %% Run simulation in loop due to memory constraints
    if simtype == 1
        [t, x, m] = tauMethod(stoich_matrix, @propensities_2state, tspan, x0, k, [], twindow/tscaleps, tscaleps);
    elseif simtype == 2
        %%% returns the timestamp, mol concentrations, and reaction (ID) that
        %%% took place at each time frame
        [t, x, x0_new, m, b, c, mol_name_new, store] = directMethod(stoich_matrix, @propensities_2state, tspan,x, x0, k, [], maxnumreactest, bond0, config0, numatoms,mol_name, store);
        
    end
    % save results and update initial conditions for next loop iteration
    tfinal(1:size(t,1),1:size(t,2)) = t; 
    xfinal(1:size(x,1), 1:size(x,2)) = x;
    mfinal(1:size(m,1),1:size(m,2)) = m;
    time = size(t,1);
    tspan = [tspan(2), tspan(2) + twindow];
    x0 = x0_new;
    bond0 = b;
    config0 = c;
    mol_name = mol_name_new;

    numt = twindow;
    while tspan(2) <= tscaleps * numframes
        if simtype == 1
            %%% calls propensities function below
            [t, x, m] = tauMethod(stoich_matrix, @propensities_2state, tspan, x0, k, [], twindow/tscaleps, tscaleps);
        elseif simtype == 2
            [t, x, x0_new, m, b, c, mol_name_new, store] = directMethod(stoich_matrix, @propensities_2state, tspan, x(end, :), x0, k, [], maxnumreactest, bond0, config0, numatoms, mol_name, store);
        end
        % save results and update initial conditions for next loop iteration
        
        tfinal(time+1:time + size(t,1),1:size(t,2)) = t; 
        xfinal(time+1:time + size(x,1), 1:size(x,2)) = x;
        mfinal(time+1:time + size(m,1),1:size(m,2)) = m;
        time = time + size(t,1);
        tspan = [tspan(2), tspan(2) + twindow];    %% shifts the time span
        x0 = x0_new;  %% final state becomes new initial conditions
        bond0 = b;
        config0 = c;
        numt = numt + twindow;  %% tracks the final time window
        mol_name = mol_name_new;
        
    end

    if tspan(1) < tscaleps * numframes    %% makes sure everything between tspan(1) and numt is computed
        tspan(2) = tscaleps * numframes;
        
        if simtype == 1
            [t, x, m] = tauMethod(stoich_matrix, @propensities_2state, tspan, x0, k, [], twindow/tscaleps, tscaleps);
        elseif simtype == 2
            [t, x, x0_new, m, b, c, mol_name_new, store] = directMethod(stoich_matrix, @propensities_2state, tspan, x(end, :), x0, k, [], maxnumreactest, bond0, config0, numatoms, mol_name, store);
        end

        % save results
        
        tfinal(time+1:time + size(t,1),1:size(t,2)) = t; 
        xfinal(time+1:time + size(x,1), 1:size(x,2)) = x;
        mfinal(time+1:time + size(m,1),1:size(m,2)) = m;
        time = time + size(t,1);
        numt = tscaleps * numframes;  %% is this not redundant with before?
        mol_name = mol_name_new;
        
    end
        
    if size(find(tfinal(2:end,1)==0))>0
        a = min(find(tfinal(2:end,1)==0));
        tfinal = tfinal(1:a-1, :);
        xfinal=xfinal(1:a-1,:);
        mfinal = mfinal(1:a-1,:);
    end
    

end

function a = propensities_2state(x, p)
    reactc = reactconc(x);  %% calls reactconc function below
    a = p.*reactc;
end

function reactc = reactconc(x)
    global ireactidx;   %% linear indices where value != 0 (column vector)
    global ireactm;     %% linear arrangement of non-zero values (column vector)
    global concexp;     %% keeps track of multiple appearances of same species

    xr = ones(size(concexp)); % default 1-0=1
    xr(ireactidx) = x(ireactm);
    xr = xr - concexp;  %% accounts for subtraction in choose function
    reactc = prod(xr, 2);    %% multiplies propensities together
    reactc(reactc < 0) = 0;



end
