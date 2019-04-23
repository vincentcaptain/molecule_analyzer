tic
close all;
clear variables;
fclose('all');

global mC1
global rpfC1
global mmC1
global newmoldict

% Reset workspace
addpath('C:\Users\Vincent\Stanford\Research\Code\Analysis\New_analysis_development\Analysis2.0\Modifying_2\matlab\');
codedir = pwd;
databasedir = 'C:\Users\Vincent\Stanford\Research\Code\Analysis\New_analysis_development\Analysis2.0\Modifying_2\matlab\input_isobutane';

%% Settings
tag = 'fullstoch'; % label for saving files
npttrain = 1; % ID of molecular dynamics simulation
b = 1; % tau = b*0.012 picoseconds
molnames = {'CH4', 'H2'}; % printed name for molecules
molfullnames = {'C1 H4 4(H-C) ', 'H2 1(H-H) '}; % name for molecules in moleculedict.txt
molsfound = length(molfullnames);
numUQ = 1; % number of Gillespie simulations
numatoms = 200;

% plot settings
plotkflag = 1;
simulflag = 1;

%% Calculated Settings
tau = [4];
% tau = [1, 4, 8, 12, 16, 20, 24, 40, 80, 120];
tau2 = tau * 0.012;
timesteps = length(tau);
error_train_iso = zeros(timesteps, molsfound);
meanmols_train_iso = zeros(timesteps, molsfound);
fluct_train_iso = zeros(timesteps, molsfound, numUQ);

for tstep = 1:timesteps
tscaleps = 0.012; % tau in picoseconds
tau_value = tau(tstep);

% parameter estimation settings
% startcalcframe = ceil(20/tscaleps) + 1; % remove (at least) first 10ps of mD for vibrational modes to equilibriate
startcalcframe = 2;
tmax = floor(0/tscaleps); % max timestep to calculate reaction rates from; default 0 = full
% tmax = ceil(100/tscaleps) + 1;

% stochastic simulation time series settings
startframe = startcalcframe - 1; % NOTE: 0-INDEXED; recall mD data start at t=1 rather than t=0

%% Get training data
datadir = [databasedir, num2str(tau(tstep))]; %% folder where data is stored
[mmD, mD, rpfC, rC, rbC] = getData(datadir); %%mD = nb of mol at each time, rpfC = nb of reaction at each time, rC = coeff for each mol for each reaction, rbC = coeff for each reactant at each time
newmoldict = getNewMolDict(datadir);
numframes = size(mD, 1);
[numreacts, nummols] = size(rC);

molids = zeros(molsfound, 1);
mmolids = zeros(molsfound, 1);
for molidx = 1:molsfound
    mmolids(molidx) = getMolidByName(molfullnames(molidx), datadir, 1);
end

%% Estimate Rates (k) for Elementary Reactions
%%% global variables to pass into Gillespie simulation
global reactset;    %% full reaction coefficients
global ireactmat;   %% finds index where non-zero and fills into front of matrix
global concexp;     %% keeps track of multiple appearances of same species
global ireactidx;   %% linear indices where value != 0 (column vector)
global ireactm;     %% linear arrangement of non-zero values (column vector)

reactset = rC;  %% numreacts by nummols, coeff of the molecule for each reaction
numreactants = full(sum(rbC, 2));  %% sum of each row in column vector numreacts by 1
maxnumreactants = max(numreactants);    %% 8 is maxsum of rows
ireactmat = zeros(numreacts, maxnumreactants);    %% numreacts by max
concexp = zeros(numreacts, maxnumreactants); % extras are 0's
cfactor = ones(numreacts, 1);
for r = 1:numreacts
    midx = 0;
    for m = find(rbC(r, :) > 0)
        for j = 1:rbC(r, m)
            midx = midx + 1;
            ireactmat(r, midx) = m;  %% pushes all indices to the front of the row
            concexp(r, midx) = j - 1;  %% counts multiple appearances
        end
        cfactor(r) = cfactor(r)*factorial(rbC(r, m));    %% Ex: 2001 + 2001 --> ..., you'll need to divide by cfactor = 2 because otherwise you will count twice the way of making this reaction happen
    end
end
ireactidx = find(ireactmat ~= 0); %% finds linear indices where value != 0 and puts into column vector
ireactm = ireactmat(ireactmat ~= 0);  %% finds all values != 0 and puts into column vector

% save mD and rpfC into global variables for elemReactGSSA_calcrates; 
% NOTE: THIS IS LEGACY CODE FOR MULTIPLE MD; MAY BE MODIFIED IN FUTURE VERSIONS
mC1 = full(mD);
rpfC1 = rpfC;
mmC1 = full(mmD);

% compute reaction rates
if tmax     %% usually false because tmax == 0
    endframe = tmax;
else
    endframe = size(rpfC, 1);
end

[k_train_iso, rarereacts,sumr,lenyesr,sumc, dk] = elemReactGSSA_calcrates(tag, 1, rC, rbC, startcalcframe, endframe, tscaleps, 1, 0, tag); %% this gets all the reaction rates.
disp('Get reaction rates');

%% Gillespie Stochastic Simulation

mDinit = zeros(nummols, 1);
[bond, config] = getBondData(datadir, numatoms ,startcalcframe);
if startframe == 0
    mDinit(molids(4)) = 56;
else
    mDinit = round(mD(startframe, :)');  %% initial values (int) in column vector
end

% run simulation
runframes = numframes - startframe; % runframes * tscaleps = total runtime
xinterp = zeros(runframes + 1, molsfound); %% creates frames by mols matrix
xinterp_avg_train_iso = zeros(runframes + 1, molsfound);  %% averages the simulations
for simulcount = 1:numUQ  %% how many Gillespie simuations we're running
    disp('Starting KMC');
   [tgssa, x, simreacts, numt, molecule_name] = elemReactGSSA_wrapper(bond, config,numatoms, mDinit, k_train_iso, 2, runframes, tscaleps); %% this actually does the advancing
   disp('Finished KMC');
   truerunframes = round(numt/tscaleps);    %% true because it's based on simulation; always <= runframes
   tscaled = tgssa/tscaleps;    %% tgssa is all timeframes concatenated into column vector
   tscaled(end) = round(tscaled(end));  %% get the last timestep
   newmolid = zeros(molsfound,1);
   for j = 1:molsfound
       for k =1:size(molecule_name, 1)
           if molfullnames(j) == molecule_name(k,1)
               newmolid(j,1) = k;  
           end
       end
       xinterp(1:truerunframes + 1, j) = interp1(tscaled, x(:, newmolid(j,1)), 0:truerunframes);    %% linear interpolated
       fluct_train_iso(tstep, j, simulcount) = sqrt(immse(xinterp(:, j), xinterp_avg_train_iso(:, j)));
   end
   xinterp_avg_train_iso = xinterp_avg_train_iso + xinterp;
   simulcount;   %% print for sanity check
end

xinterp_avg_train_iso = xinterp_avg_train_iso/numUQ;
trange = 1 + (0:truerunframes);   %% 1-indexed row vector of runframes
timerange = (startframe + trange - 1) * tscaleps; %% row vector of timestamps
timerangeidx = startframe + (0:runframes);
timerange = timerangeidx * tscaleps;

% save('all_data.mat','-v7.3');

% Plot Simulation Results
if simulflag
    for merr = 1:molsfound
        figure(merr); hold all;
        plot(timerange', xinterp_avg_train_iso(:, merr), 'linewidth', 3.5);
    end
end

nnzk = sum(k_train_iso~=0);

end

%% Plot the MD simulations
for merr = 1:molsfound
    figure(merr);
    molidx = mmolids(merr);
    plot(timerange', mmC1(timerangeidx, molidx), 'k', 'linewidth', 3);
    h = legend('KMC', 'MD', 'fontsize', 24, 'Location', 'eastoutside');
    set(gca, 'fontsize', 30);
    set(gca, 'ticklength', [0.03 0.03]);
    set(gca, 'linewidth', 3.5);
    xlabel('Time (ps)');
    ylabel('Number of molecules');
    title(molnames(merr));
end

%% Save figures to file
% for fig = 1:molsfound+2
%     print(figure(fig), ['Figure_train_iso_', num2str(fig)], '-depsc', '-r600');
%     print(figure(fig), ['Figure_train_iso_', num2str(fig)], '-dpdf', '-r600');
%     print(figure(fig), ['Figure_train_iso_', num2str(fig)], '-dpng', '-r600');
% end
% save('fluct_train_iso.mat', 'fluct_train_iso');
% save('error_train_iso.mat', 'error_train_iso', 'meanmols_train_iso')
save('test_2.dat');
toc

