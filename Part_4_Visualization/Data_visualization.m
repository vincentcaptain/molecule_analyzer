close all;
clear variables;
fclose('all');

K = load('all_data.mat', '-mat');

%%
datadir = '~/Downloads/input_isobutane4';
molsfound= K.molsfound;
timerange = K.timerange;
xinterp_avg_train_iso = K.xinterp_avg_train_iso;
mmolids = K.mmolids;
mmC1 = K.mmC1;
timerangeidx = K.timerangeidx;
molnames = K.molnames;
molecule_name = K.molecule_name;
x = K.x;
k_train_iso = K.k_train_iso;
runframes = K.runframes;
tscaled = K.tscaled;
truerunframes = K.truerunframes;
numUQ = K.numUQ;

%%
% 
% xinterp_avg_train_iso = K.xinterp_avg_train_iso;
% figure(1); hold all;
% plot(timerange', xinterp_avg_train_iso(:, 1), 'linewidth', 3.5);



%% Plot the molecules of interest
close all;


molnames = {'CH4', 'H2', 'C2H6'}; % printed name for molecules
molfullnames = {'C1 H4 4(H-C) ', 'H2 1(H-H) ', 'C2 H6 1(C-C) 6(H-C) '}; % name for molecules in moleculedict.txt
molsfound = length(molfullnames);

mmolids = zeros(molsfound, 1);
for molidx = 1:molsfound
    mmolids(molidx) = getMolidByName(molfullnames(molidx), datadir, 1);
end


xinterp = zeros(runframes + 1, molsfound); %% creates frames by mols matrix
xinterp_avg_train_iso = zeros(runframes + 1, molsfound);  %% averages the simulations
newmolid = zeros(molsfound,1);
for j = 1:molsfound
   for k =1:size(molecule_name, 1)
       if molfullnames(j) == molecule_name(k,1)
           newmolid(j,1) = k;  
       end
   end
   xinterp(1:truerunframes + 1, j) = interp1(tscaled, x(:, newmolid(j,1)), 0:truerunframes);    %% linear interpolated
end
xinterp_avg_train_iso = xinterp_avg_train_iso + xinterp;
xinterp_avg_train_iso = xinterp_avg_train_iso/numUQ;


% Plot Simulation Results
%     if tau(tstep) == 1 || tau(tstep) == 8 || tau(tstep) == 16
        for merr = 1:molsfound
            figure(merr); hold all;
            plot(timerange', xinterp_avg_train_iso(:, merr), 'linewidth', 3.5);
        end
%     end

% nnzk = sum(k_train_iso~=0);
% % save(['xinterp_train_iso_', num2str(tau(tstep)), '.mat'], 'timerangeidx', 'timerange', 'xinterp_avg_train_iso', 'k_train_iso')
% 
% % Find rms error
% for merr = 1:molsfound
%     molidx = molids(merr);
%     error_train_iso(tstep, merr) = sqrt(immse(xinterp_avg_train_iso(:, merr), mC1(timerangeidx, molidx)));
%     meanmols_train_iso(tstep, merr) = mean(xinterp_avg_train_iso(:, merr));
% end

% Plot the MD simulations
for merr = 1:molsfound
    figure(merr);
    molidx = mmolids(merr);
    plot(timerange', mmC1(timerangeidx, molidx), 'k', 'linewidth', 3);
    h = legend('KMC', 'MD', 'fontsize', 24, 'Location', 'eastoutside');
%     P = get(h, 'Position');
%     hAxes = axes('Parent', h.Parent, 'Units', h.Units, 'Position', [P(1)-0.005 P(2) P(3) P(4)], 'XTick', [], 'YTick', [], ...
%         'Color', 'none', 'YColor', 'none', 'XColor', 'none', 'HandleVisibility', 'off', 'HitTest', 'off');
%     hTitle = title(hAxes, '\tau values');
%     set(hTitle, 'fontsize', 28);
    set(gca, 'fontsize', 30);
    set(gca, 'ticklength', [0.03 0.03]);
    set(gca, 'linewidth', 3.5);
    xlabel('Time (ps)');
    ylabel('Number of molecules');
    title(molnames(merr));
end
%% Plot carbon clusters for MD
[biggest_CC, num_CC, numofCinCC] = getBiggestCarbonCluster(mmC1, datadir, 1, []);

figure(10)
plot(timerange',biggest_CC(timerangeidx), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Size of biggest carbon cluster');
title('MD');

figure(11)
plot(timerange',num_CC(timerangeidx), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of carbon cluster');
title('MD');

figure(12)
plot(timerange',numofCinCC(timerangeidx), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of carbons in carbon clusters');
title('MD');

%% Plot carbon clusters for KMC
[KMC_biggest_CC, KMC_num_CC, KMC_numofCinCC] = getBiggestCarbonCluster(x, datadir, 0, molecule_name);

figure(13)
plot(KMC_biggest_CC(:), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Size of biggest carbon cluster');
title('KMC');

figure(14)
plot(KMC_num_CC(:), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of carbon clusters');
title('KMC');

figure(15)
plot(KMC_numofCinCC(:), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of carbon in carbon clusters');
title('KMC');

%% Plot reaction rates
figure(100);
h = stem(k_train_iso, ':diamondr', 'fill');
hbase = get(h, 'Baseline');
set(hbase, 'Visible', 'off');
xlabel('Elementary Reaction Index', 'FontSize', 16);
ylabel('k (conc/picosecond)', 'FontSize', 16);
title('Reaction Rate Coefficients', 'FontSize', 20);
set(gca, 'FontSize', 14);

%% Get the molecules for one frame
% frame = 9866;
% num_atoms = size(x,2);
% molhist = find(x(frame,:)>0);
% size(molhist)
% for i = 1:size(molhist,2)
%     [molecule_name(molhist(i)), x(frame, molhist(i))]
% end

%%

[num_cycles, numofCincycles, numofHincycles, diff] = getCycles(mmC1, datadir, 1, molecule_name);


figure(16)
plot(timerange',num_cycles(timerangeidx), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of cycles');
title('MD');

% figure(17)
% plot(timerange',numofCincycles(timerangeidx), 'k', 'linewidth', 3);
% xlabel('Time (ps)');
% ylabel('Number of carbon in cycles');
% title('MD');
% 
% figure(18)
% plot(timerange',numofHincycles(timerangeidx), 'k', 'linewidth', 3);
% xlabel('Time (ps)');
% ylabel('Number of hydrogens in cycles');
% title('MD');

figure(23)
plot(timerange',diff(timerangeidx), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Diff');
title('MD');

[KMC_num_cycles, KMC_numofCincycles, KMC_numofHincycles, KMC_diff] = getCycles(x, datadir, 0, molecule_name);

figure(19)
plot(KMC_num_cycles(:), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Number of cycles');
title('KMC');
% 
% figure(20)
% plot(KMC_numofCincycles(:), 'k', 'linewidth', 3);
% xlabel('Time (ps)');
% ylabel('Number of carbon in cycles');
% title('KMC');
% 
% figure(21)
% plot(KMC_numofHincycles(:), 'k', 'linewidth', 3);
% xlabel('Time (ps)');
% ylabel('Number of hydrogen in cycles');
% title('KMC');

figure(22)
plot(KMC_diff(:), 'k', 'linewidth', 3);
xlabel('Time (ps)');
ylabel('Diff');
title('KMC');

%%

