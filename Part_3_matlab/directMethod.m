function [t, x, x0_new, m, b, c, mol_name_final, store] = directMethod(stoich_matrix, propensity_fcn, tspan,x_mol,  x0,...
                                  rate_params, output_fcn, MAX_OUTPUT_LENGTH, bond0, config0, numatoms,mol_name, store)
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
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, size(x0,2)+10);
MU = zeros(MAX_OUTPUT_LENGTH, 1);
T(1) = tspan(1);    %% initial time value
B = bond0;
C = config0;
rxn_count = 1;

%% MAIN LOOP
while T(rxn_count) < tspan(2)
    % Calculate reaction propensities
    a  = propensity_fcn(x0, rate_params);   %% i-th reaction
    % Compute tau and mu using random variates
    a0 = sum(a);
    r = rand(1, 2);  
    tau = (1/a0) * log(1 / r(1));   
    mu  = find((cumsum(a) >= r(2) * a0), 1, 'first');
    while isPossible(mu, C, B, store) ~= 1
        r = rand(1, 2);  
        tau = (1/a0) * log(1 / r(1));   
        mu  = find((cumsum(a) >= r(2) * a0), 1, 'first');
    end
    
    % Update time and carry out reaction mu
    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count, :);
        b = B;
        c = C;
        disp('Be careful youre here');
        warning('SSA:ExceededCapacity',...
                ['Number of reaction events exceeded the number pre-allocated (',num2str(MAX_OUTPUT_LENGTH),'. Simulation terminated prematurely.']);
        return;
    end
    
    [x_new, x0, bond_new, config_new, mol_name_new, store] = makeReaction(B, C,mu, mol_name, x_mol, x0, store);
    mol_name = mol_name_new;
    X(rxn_count, 1:size(x_new,2)) = x_new; 
    MU(rxn_count)  = mu;
    B = bond_new;
    C = config_new;
    
    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count, :)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count, :);
            m = MU(1:rxn_count);
            b = B;
            c = C;
            x0_new = x0;
            disp('Ayayaya');
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end
    x_mol = x_new;
    T(rxn_count + 1) = T(rxn_count) + tau;
    rxn_count = rxn_count + 1;
end  

X = X(:, 1:size(x_new,2));


% Record output
t = T(1:rxn_count-1);
x = X(1:rxn_count-1, :);
m = MU(1:rxn_count-1);
b = B;
c = C;
x0_new = x0;
mol_name_final = mol_name;

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
