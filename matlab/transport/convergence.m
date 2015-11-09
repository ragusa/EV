close all; clc;

%% User Options
%--------------------------------------------------------------------------
% mesh options
%--------------------------------------------------------------------------
nel_init = 50; % number of elements in first cycle
n_cycle  = 1;  % number of refinement cycles
%--------------------------------------------------------------------------

%% Setup

% mesh sizes and L2 errors for each cycle
h        = zeros(n_cycle,1);
uL_err   = zeros(n_cycle,1);
uH_err   = zeros(n_cycle,1);
uFCT_err = zeros(n_cycle,1);

%% Perform runs
for cycle = 1:n_cycle
    nel = nel_init * 2^(cycle-1); % number of elements
    n_dof = nel + 1;              % number of dofs
    h(cycle) = 1/nel;             % mesh size normalized by length

    fprintf('Cycle %i: n_dof = %i\n',cycle,n_dof);
    
    % perform runs
    
    % compute error for runs
end

%% Convergence Table

fprintf('\nLow-order Convergence:\n')
for cycle = 1:n_cycle
    if (cycle == 1)
        fprintf('Cycle %i: L2 error = %e\n',cycle,uL_err(cycle));
    else
        rate = log(uL_err(cycle)/uL_err(cycle-1))/log(h(cycle)/h(cycle-1));
        fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uL_err(cycle),rate);
    end
end
fprintf('\nHigh-order Convergence:\n')
for cycle = 1:n_cycle
    if (cycle == 1)
        fprintf('Cycle %i: L2 error = %e\n',cycle,uH_err(cycle));
    else
        rate = log(uH_err(cycle)/uH_err(cycle-1))/log(h(cycle)/h(cycle-1));
        fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uH_err(cycle),rate);
    end
end
fprintf('\nFCT Convergence:\n')
for cycle = 1:n_cycle
    if (cycle == 1)
        fprintf('Cycle %i: L2 error = %e\n',cycle,uFCT_err(cycle));
    else
        rate = log(uFCT_err(cycle)/uFCT_err(cycle-1))/log(h(cycle)/h(cycle-1));
        fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uFCT_err(cycle),rate);
    end
end

%% Convergence Plot