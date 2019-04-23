function [t, x, x0_new, m, b, c, mol_name_final] = directMethod(stoich_matrix, propensity_fcn, tspan, x0,...
                                  rate_params, output_fcn, MAX_OUTPUT_LENGTH, bond0, config0, numatoms,mol_name)
%DIRECTMETHOD Implementation of the Direct Method variant of the Gillespie algorithm
%   Based on: Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%   Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%
%   Usage:
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0 )
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0, rate_params )
%       [t, x] = directMethod( stoich_matrix, propensity_fcn, tspan, x0, rate_params, output_fcn )
%
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)    
%
%   Required:
%       tspan:          Initial and final times, [t_initial, t_final].
%       x0:             Initial species values, [S1_0, S2_0, ... ].
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%       prop_fcn:       Function handle to function that calculates
%                       reaction propensities.
%                       Target function should be of the form
%                           ac = f( xc, rate_params ),
%                       where xc is the current state [S1, S2, ...], and
%                       rate_params is the user-defined rate parameters.
%                       The function should return vector ac of 
%                       propensities (Nreactions x 1) in the same order as
%                       the reactions given in stoich_matrix.
%
%	Optional:
%		output_fcn:	    Handle to user-defined function with the form
%							status = f( tc, xc )
%						The output_fcn is called at each time step in the 
%						simulation and is passed the current time and state.
%						It can be used to locate events or monitor progress. 
%						If it returns 0, the simulation terminates.
%
%   Author:  Nezar Abdennur, 2012 <nabdennur@gmail.com>
%   Created: 2012-01-19
%   Dynamical Systems Biology Laboratory, University of Ottawa
%   www.sysbiolab.uottawa.ca
%
%   See also FIRSTREACTIONMETHOD
%
%   Modified: 2014-08-12
%   Modified By: Qian Yang <qianyang@stanford.edu>
%   Now also outputs reactions in vector m.

if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end
if ~exist('rate_params', 'var')
    rate_params = [];
end
    
%% Initialize
%num_rxns = size(stoich_matrix, 1);
% num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, size(x0,2)+10);
MU = zeros(MAX_OUTPUT_LENGTH, 1);
B = zeros(1 ,numatoms, 10);
C = zeros(1, numatoms);
T(1) = tspan(1);    %% initial time value
B(1,:, :) = bond0;
C(1,:) = config0;
% if size(mol_name,1) == 1
%     X = zeros(1,1);
% else
%     [X, mol_name_new] = getMoleculePerFrame(B(1,:,:),C(1,:), mol_name);        %% initial concentration value
% end
rxn_count = 1;
% mol_name = mol_name_new;

%% MAIN LOOP
while T(rxn_count) < tspan(2)
    % Calculate reaction propensities
%     a  = propensity_fcn(X(rxn_count, :), rate_params);   %% i-th reaction
    a  = propensity_fcn(x0, rate_params);   %% i-th reaction
    % Compute tau and mu using random variates
    a0 = sum(a);
    r = rand(1, 2);  
    tau = (1/a0) * log(1 / r(1));   
    mu  = find((cumsum(a) >= r(2) * a0), 1, 'first');
    while isPossible(mu, C, B) ~= 1
        r = rand(1, 2);  
        tau = (1/a0) * log(1 / r(1));   
        mu  = find((cumsum(a) >= r(2) * a0), 1, 'first');
    end
    % alternatively...
    %mu=1; s=a(1); r0=r(2)*a0;
    %while s < r0
    %   mu = mu + 1;
    %   s = s + a(mu);
    %end
    
    % Update time and carry out reaction mu
    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count, :);
        b = B(1, :, :);
        c = C(1, :);
        disp('Be careful youre here');
        warning('SSA:ExceededCapacity',...
                ['Number of reaction events exceeded the number pre-allocated (',num2str(MAX_OUTPUT_LENGTH),'. Simulation terminated prematurely.']);
        return;
    end
    
%     T(rxn_count + 1) = T(rxn_count) + tau;
    
    
%     if mod(rxn_count,100000)==0
%         test = T(rxn_count+1)
%     end
    
%     if mu < 1
%         mu
%     end
%     if length(X(rxn_count, :))~=length(stoich_matrix(mu, :))
%         mu
%     end
%     if isempty(stoich_matrix(mu, :))
%         mu
%     end

    % HERE BIG THING TO MODIFY !!!!!!!!!!
%     if rxn_count >1
%         disp('here')
%         find(C(rxn_count,:) - C(rxn_count-1,:) ~= 0)
%         find(B(rxn_count,:,:) - B(rxn_count-1,:, :) ~= 0)
%     end
    [x_new, x0, bond_new, config_new, mol_name_new] = makeReaction(B(1,:,:), C(1,:),mu, mol_name);
    mol_name = mol_name_new;
    
    X(rxn_count, 1:size(x_new,2)) = x_new; 
    
%     X2 = zeros(rxn_count,max(size(x_new,2), size(X,2)));
%     if rxn_count>1
%         X2(1:rxn_count-1,1:size(X,2)) = X;
%     end
%     X2(rxn_count, 1:size(x_new,2)) = x_new;
%     X = X2;

%     X(rxn_count + 1, :) =  x_new;
    MU(rxn_count)  = mu;
%     B(rxn_count +1, :, :) = bond_new;
%     C(rxn_count + 1, :) = config_new;
%     B(2, :,:) = B(1,:,:);
    B(1,:,:) = bond_new;
%     C(2,:) = C(1,:);
    C(1,:) = config_new;
    
    
%     rxn_count = rxn_count + 1;
    
    
%     if sum(X(rxn_count, :) < 0)
%         mu
%     end

    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count, :)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count, :);
            m = MU(1:rxn_count);
            b = B(1, :, :);
            c = C(1, :);
            x0_new = x0;
            disp('Ayayaya');
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end
    
    T(rxn_count + 1) = T(rxn_count) + tau;
    rxn_count = rxn_count + 1;
end  
% disp(['rxn_count: ',num2str(rxn_count)]);

X = X(:, 1:size(x_new,2));


% Record output
t = T(1:rxn_count-1);
x = X(1:rxn_count-1, :);
m = MU(1:rxn_count-1);
b = B(1, :, :);
c = C(1, :);
x0_new = x0;
mol_name_final = mol_name;

% 
% if t(end) > tspan(2)
%     t(end) = tspan(2);
%     x(end, :) = X(rxn_count - 1, :);
%     m(end) = 0;
%     b = B(2, :, :);
%     c = C(2, :);
%     
%     [x0_new, ] = getNumOfMol(b, c, datadir);
%     
% %     find(test-x0 ~=0)
% %     test(find(test-x0 ~=0))
% %     x0(find(test-x0 ~=0))
% 
% end    


end

function stoicherror = checkstoich(molcount)

    global molatomdict;
    
    s = molcount * molatomdict;
    if sum(molcount < 0) ~= 0
        stoicherror = 3;
    elseif s(1) ~= 216
        stoicherror = 1;
    elseif s(2) ~= 216 * 4
        stoicherror = 2;
    else
        stoicherror = 0;
    end
end
