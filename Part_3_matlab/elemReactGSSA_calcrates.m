function [k,rarereacts,sumr,lenyesr,sumc, dk] = elemReactGSSA_calcrates(kname,setlist,mDrb,mDrpos,startframe,endframe,tscaleps,rarereactlim,computeksave,computekname)
    %%% tag is kname and computekname, 1 is setlist, rC is mDrb, rbC is
    %%% mDrpos, rarereactlim and computeksave are both 0.

    global mC1; % mC2 mC3 mC4 mC5 mC6 mC7 mC8 mC9 mC10 mC11 mC12
    global rpfC1; % rpfC2 rpfC3 rpfC4 rpfC5 rpfC6 rpfC7 rpfC8 rpfC9 prfC10 rpfC11 rpfC12
    global ireactmat;   %% finds index where non-zero and fills into front of matrix
    global concexp;     %% keeps track of multiple appearances of same species
  
    k = zeros(size(mDrb, 1), 1);  %% number of reactions
    backt = 1;
    numtimesteps = endframe - startframe + 1;   %% probably endframe == numframes
    
    dk = zeros(size(k));
    sumr = zeros(size(k));  %% same sizes as k
    sumc = zeros(size(k));
    lenyesr = zeros(size(k));
    
    for setidx = setlist  %% currently setlist = 1
        
        for r = 1:size(eval(['rpfC', num2str(setidx)]), 2)  %% r = 1:numreactions
            % dk/dt at each time point is mDrpf(t,r)
            dkdt = eval(['rpfC', num2str(setidx), '(startframe:endframe, r)']);    %% evaluates one reaction across all time frames
            if sum(dkdt) == 0 %% that reaction never occurs(?)
                continue;
            end
            
            % order of mols in react are in exponents of mDrpos
            ireact = find(mDrpos(r, :) > 0); %% gets (column) ID of reactants in row vector
            concexpireact = mDrpos(r, ireact);   %% now returns those coefficients from the IDs
            concexpireactmat = ones(length(1:numtimesteps - backt), 1)*concexpireact;  %% duplicates coefficients for all timestep-rows

            % reactions can only occur if all reactants enough concentration;
            % do all calculations in matrix form for speedup
            %%% compares the molecule amounts with stoich in reaction
            reactnotready = sum(eval(['mC', num2str(setidx), '(startframe:endframe-backt,ireact)']) < concexpireactmat, 2);
%             reactnotready = sum(mD(1:end-1,ireact)<concexpireactmat,2);
            reactready = find(reactnotready == 0);    %% if sum to 0 then all are >= threshold
            yesreact = reactready + backt; % reaction timestep is timestep after concentration available. why??
            
            xr = ones(length(yesreact), sum(concexpireact)); %% sum covers all possible values in row
            %%% scales back to original mD table, 
            xr = eval(['mC', num2str(setidx), '(reactready+startframe-1,ireactmat(r,1:sum(concexpireact)))']); %% when the reaction is ready, store in xr the number of molecules of each reactant
            xr = xr - ones(length(yesreact), 1) * concexp(r, 1:sum(concexpireact)); %% accounts for subtraction in choose function, if you have twice the same reactant, the first time you have n choices but the 2nd time only n-1 and etc...
            concreact = prod(xr, 2); %% multiplies all values in row together (propensity function), how many ways of doing a reaction
            concreact(concreact<0) = 0; %% take all negative numbers and set to 0, maybe if you said reaction is ready but 
%             concreact = prod(eval(['mC',num2str(setidx),'(reactready+startframe-1,ireact)']).^concexpireactmat(1:length(yesreact),:),2);
            

            dk(r) = dk(r) + sum(dkdt);
            sumr(r) = sumr(r) + sum(dkdt(yesreact));    %% records total number of times reaction appears
            lenyesr(r) = lenyesr(r) + length(yesreact); %% records how many instances of good reactions
            sumc(r) = sumc(r) + sum(concreact) + dk(r) - sumr(r);         %% records sum of combinations (propensity functions)
        end
        
        setidx;   %% print for a sanity check
    end

    % defined by MLE of Poisson rv
    % There is a modification here compared with Qian's code and
    % calculations: in the local case, it is possible to have
    % several reactions happening at one step. However, the
    % reactants of the 2nd reaction (for example) may not be
    % present at the beginning of the timestep, so yesreact = 0. In
    % order to deal with that, I consider that the reaction was ready if it
    % happened.
    k = (dk./(sumc))/(tscaleps);  %% equation 4 in Qian's paper
    k(sumc == 0) = 0; %% where concentrations are 0, also set k to 0. Only for combining data sets
    % remove rare events
    k(dk < rarereactlim) = 0;   %% Set to zero, but naive threshold for minimum occurrences
    rarereacts = find(sumr < rarereactlim);   %% indices that are rare
%         if length(yesreact) < rarereactlim
%             disp(['Likely overestimate: ',num2str(k(r))]);
%             k(r) = 0;
%             numrarereact = numrarereact + 1;
%         end
    if sum(k > 1/tscaleps)
       disp('Error: some k-values are too large');
       return;
    end
    if isnan(k)
        disp('Error: contains NaN k');
        return;
    end
    
    if computeksave
        save([kname, '_', computekname], 'k', '-ascii');
    end
end
