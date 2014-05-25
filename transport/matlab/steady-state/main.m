close all; clear; clc;
% Solves a steady-state transport equation:
%   -(K  )*u = b (high-order system)
%   -(K+D)*u = b (low-order system)

nel_init = 10; % number of elements in first cycle
n_cycle  = 4;  % number of refinement cycles

problemID = 1; % 1 = smooth exponential decay
               % 2 = sharp exponential decay; heterogenous sigma

limiting_option = 2; % 0 = set limiting coefficients to 0 (no correction)
                     % 1 = set limiting coefficients to 1 (full correction)
                     % 2 = compute limiting coefficients normally
max_iter = 10; % maximum number of iterations for implicit FCT
%==========================================================================

% determine problem parameters from ID
switch problemID
    case 1
        len=10;         % length of domain
        omega=1.0;      % omega
        sigma = @(x) 1; % sigma
        src=0;          % source
        inc=1;          % incoming flux
    case 2
        len=10;
        omega=1.0;
        sigma = @hetero_sigma;
        src=0;
        inc=1;
    otherwise
        error('Invalid problem ID');
end

h = zeros(n_cycle,1); % mesh size for each cycle

% L2 errors for each cycle
uL_err   = zeros(n_cycle,1);
uH_err   = zeros(n_cycle,1);
uFCT_err = zeros(n_cycle,1);

% loop over refinement cycles
for cycle = 1:n_cycle
    nel = nel_init * 2^(cycle-1); % number of elements
    n_dof = nel + 1;              % number of dofs
    h(cycle) = len/nel;           % mesh size
    
    % build matrices
    [MC,ML,K,D,b] = build_matrices(len,nel,omega,sigma,src);
    ML = diag(ML);
    AH = -(K  );    % high-order steady-state system matrix
    AL = -(K+D);    % low-order steady-state system matrix
    
    % compute modified high-order sytem matrix for Dirichlet BC
    AH_mod = AH;
    AH_mod(1,:)=0; AH_mod(1,1)=1;
    % compute modified low-order sytem matrix for Dirichlet BC
    AL_mod = AL;
    AL_mod(1,:)=0; AL_mod(1,1)=1;
    % compute modified rhs
    b_mod = b;
    b_mod(1) = inc;
    
    % low-order solve
    %======================================================================
    fprintf('Computing low-order solution...\n');
    uL = AL_mod \ b_mod;
    
    % high-order solve
    %======================================================================
    fprintf('Computing high-order solution...\n');
    uH = AH_mod \ b_mod;
    
    % FCT solve
    %======================================================================
    fprintf('Computing FCT solution...\n');
    uFCT = uL;
    % compute the upper and lower bound for max principle
    [W_max,W_min] = compute_max_principle_bounds(uL,AL,b);
    for iter = 1:max_iter
        % FCT solve
        %----------------
        % compute flux correction matrix
        F = compute_kuzmin_flux_correction_matrix(uH,D);
        % compute limiting coefficients
        switch limiting_option
            case 0 % no correction
                limiter = zeros(n_dof,n_dof);
            case 1 % full correction (no limiting)
                limiter = ones(n_dof,n_dof);
            case 2 % normal limiting
                limiter = compute_limiting_coefficients_kuzmin(F,uFCT,W_max,W_min,AL,b);
            otherwise
                error('Invalid limiting option');
        end
        % compute correction rhs
        rhs = b + sum((limiter.*F),2);
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        uFCT = AL_mod \ rhs;
        % check that max principle is satisfied. If not, continue iteration
        satisfied_max_principle = check_max_principle(uFCT,W_max,W_min);
        if (satisfied_max_principle)
            fprintf('Satisfied max principle at iteration %i\n',iter);
            break;
        else
            fprintf('Did not satisfy max principle at iteration %i\n',iter);
        end
    end
    
    % compute error
    %======================================================================
    x = linspace(0,len,nel+1);
    u_exact = exact_solution(problemID,x,sigma,src,inc,omega);
    uL_err(cycle)   = norm(uL   - u_exact',2);
    uH_err(cycle)   = norm(uH   - u_exact',2);
    uFCT_err(cycle) = norm(uFCT - u_exact',2);
end

% compute convergence rates
%==========================================================================
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

% plot
%==========================================================================
hold all;
% exact solution
x_exact = linspace(0,len,1000);
u_exact = exact_solution(problemID,x_exact,sigma,src,inc,omega);
plot(x_exact,u_exact,'k-');
% numerical solutions
plot(x,uL,'r-s');
plot(x,uH,'b-+');
plot(x,uFCT,'g-x');
plot(x,W_min,'m:');
plot(x,W_max,'m:');
% plot legend
legend_entries = char('Exact');
legend_entries = char(legend_entries,'Low-order');
legend_entries = char(legend_entries,'High-order');
switch limiting_option
    case 0 % no correction
        limiter_string = 'no correction';
    case 1 % full correction (no limiting)
        limiter_string = 'not limited';
    case 2 % normal limiting
        limiter_string = 'limited';
    otherwise
        error('Invalid limiting option');
end
FCT_legend_string = ['FCT, ',limiter_string];
legend_entries = char(legend_entries,FCT_legend_string);
legend_entries = char(legend_entries,'FCT Min bound');
legend_entries = char(legend_entries,'FCT Max bound');
legend(legend_entries);