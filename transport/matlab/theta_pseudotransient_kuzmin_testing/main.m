close all; clc;
%% Introduction

% Solves a transient transport equation:
%   MC*du/dt-(K  )*u = b (high-order system)
%   ML*du/dt-(K+D)*u = b (low-order system)
% using the theta-scheme for time discretization

%% User Options

nel_init = 50; % number of elements in first cycle
n_cycle  = 1;  % number of refinement cycles

problemID = 1; % 1 = smooth exponential decay
               % 2 = sharp exponential decay; heterogenous sigma

theta = 0.0; % 0:   explicit euler
             % 1:   implicit euler
             % 1/2: crank nicholson
ss_tol = 1.0e-6; % steady-state tolerance
n_max = 200;      % maximum number of time steps
CFL = 0.8;       % CFL number

limiting_option = 2; % 0 = set limiting coefficients to 0 (no correction)
                     % 1 = set limiting coefficients to 1 (full correction)
                     % 2 = compute limiting coefficients normally

save_convergence_data = 0; % option to save convergence data (0=no,1=yes)

%% Setup

% determine problem parameters from ID
switch problemID
    case 1
        len=10;         % length of domain
        omega=1.0;      % omega
        sigma = @(x) 1; % sigma
        src=0;          % source
        inc=1;          % incoming flux
        speed=1.0;      % transport speed
    case 2
        len=10;
        omega=1.0;
        sigma = @hetero_sigma;
        src=0;
        inc=1;
        speed=1.0;
    case 3
        len=10;
        omega=1.0;
        sigma = @(x) 1;
        src=0;
        inc=1;
        speed=1.0;
    otherwise
        error('Invalid problem ID');
end

h = zeros(n_cycle,1); % mesh size for each cycle

% L2 errors for each cycle
uL_err   = zeros(n_cycle,1);
uH_err   = zeros(n_cycle,1);
uFCT_err = zeros(n_cycle,1);

%% Solution

% loop over refinement cycles
for cycle = 1:n_cycle
    nel = nel_init * 2^(cycle-1); % number of elements
    n_dof = nel + 1;              % number of dofs
    h(cycle) = len/nel;           % mesh size
    fprintf('Cycle %i: n_dof = %i\n',cycle,n_dof);
    
    % build matrices
    [MC,ML,K,D,b] = build_matrices(len,nel,omega,sigma,src,speed,inc);
    ML = diag(ML);
    AL = -(K+D);    % low-order steady-state system matrix
    AH = -K;        % high-order steady-state system matrix
    
    % compute dt using CFL condition for explicit Euler, even if not
    % using explicit Euler, just for consistency in dt for comparison
    dt = CFL*ML(1,1)/AL(1,1);
    for i = 2:n_dof
        dt = min(dt, CFL*ML(i,i)/AL(i,i));
    end
    
    % define low-order and high-order transient sytem matrices to be inverted
    ALtr = ML + theta*dt*AL;
    AHtr = MC + theta*dt*AH;
    % compute modified transient sytem matrix for Dirichlet BC
    ALtr_mod = ALtr;
%     ALtr_mod(1,:)=0; ALtr_mod(1,1)=1;
    % compute modified transient sytem matrix for Dirichlet BC
    AHtr_mod = AHtr;
%     AHtr_mod(1,:)=0; AHtr_mod(1,1)=1;
    
    % initial conditions for pseudotransient (equal to exact steady-state solution)
    x = linspace(0,len,nel+1);
    u0 = zeros(n_dof,1);
    
    % low-order solve
    %======================================================================
    fprintf('\tComputing low-order solution...\n');
    u_old = u0;
    for n = 1:n_max
        fprintf('\t\tTime step %i:',n);
        % compute rhs
        rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b;
        % modify rhs for Dirichlet BC
%         rhs(1) = inc;
        % solve modified system
        uL = ALtr_mod \ rhs;
        % test steady-state convergence
        ss_err = norm(uL-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        converged = ss_err < ss_tol;
        if (converged)
            fprintf('\t\tReached steady-state at time step %i\n',n);
            break;
        end
        % reset u_old
        u_old = uL;
    end
    
    % high-order solve
    %======================================================================
    fprintf('\tComputing high-order solution...\n');
    u_old = u0;
    for n = 1:n_max
        fprintf('\t\tTime step %i:',n);
        % compute rhs
        rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
        % modify rhs for Dirichlet BC
%         rhs(1) = inc;
        % solve modified system
        uH = AHtr_mod \ rhs;
        % test steady-state convergence
        ss_err = norm(uH-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        converged = ss_err < ss_tol;
        if (converged)
            fprintf('\t\tReached steady-state at time step %i\n',n);
            break;
        end
        % reset u_old
        u_old = uH;
    end
    
    % FCT solve
    %======================================================================
    fprintf('\tComputing FCT solution...\n');
    u_old = u0;
    uFCT  = u0;
    for n = 1:n_max
        fprintf('\t\tTime step %i:',n);
        
        % compute auxiliary solution
        %----------------
        % compute rhs
        rhs = (ML - (1-theta)*dt*AL)*u_old + (1-theta)*dt*b;
        % modify rhs for Dirichlet BC
%         rhs(1) = inc;
        % solve for aux. solution
        ML_mod = ML; ML_mod(1) = 1.0;
        u_aux = ML_mod \ rhs;
    
        % high-order solve
        %----------------
        % compute rhs
        rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
        % modify rhs for Dirichlet BC
%         rhs(1) = inc;
        % solve modified system
        uH_FCT_loop = AHtr_mod \ rhs;
        
        % FCT solve
        %----------------
        % compute flux correction matrix
        F = flux_correction_matrix(u_old,uH_FCT_loop,dt,D,MC,theta);
        % cancel antidiffusive fluxes down the gradient of u_aux
        F = cancel_down_gradient(F,u_aux);
        % compute limiting coefficients
        switch limiting_option
            case 0 % no correction
                Flim = 0*F;
            case 1 % full correction (no limiting)
                Flim = F;
            case 2 % normal limiting
                Flim = limit_fluxes(F,u_aux,ML,AL,b,dt,theta);
            otherwise
                error('Invalid limiting option');
        end
        % enforce LED
        Flim = enforce_LED(Flim,u_aux);
        % compute correction rhs
        rhs = ML*u_aux + theta*dt*b + sum(Flim,2);
        % modify rhs for Dirichlet BC
%         rhs(1) = inc;
        % solve for FCT solution
        uFCT = ALtr_mod \ rhs;
        
        % check maximum principle
        satisfied_max_principle = check_max_principle(u_aux,uFCT);
        if (~satisfied_max_principle)
            error('Did not satisfy maximum principle\n');
        end

        % plot
        plot(x,uFCT);
        pause(0.1);
        
        % test steady-state convergence
        ss_err = norm(uFCT-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        converged = ss_err < ss_tol;
        if (converged)
            fprintf('\t\tReached steady-state at time step %i\n',n);
            break;
        end
        % reset u_old
        u_old = uFCT;
    end
    
    % compute error
    %======================================================================
    uL_err(cycle)   = compute_error(uL,x);
    uH_err(cycle)   = compute_error(uH,x);
    uFCT_err(cycle) = compute_error(uFCT,x);
end

%% Convergence Rates

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

if (save_convergence_data)
    % save convergence data to file
    dlmwrite('output/convergence.csv',[h,uL_err,uH_err,uFCT_err]);
end

%% Plot

% hold all;
% % exact solution
% x_exact = linspace(0,len,1000);
% u_exact = exact_solution(problemID,x_exact,sigma,src,inc,omega,speed);
% plot(x_exact,u_exact,'k-');
% % numerical solutions
% plot(x,uL,'r-s');
% plot(x,uH,'b-+');
% plot(x,uFCT,'g-x');
% % plot legend
% legend_entries = char('Exact');
% legend_entries = char(legend_entries,'Low-order');
% legend_entries = char(legend_entries,'High-order');
% switch limiting_option
%     case 0 % no correction
%         limiter_string = 'no correction';
%     case 1 % full correction (no limiting)
%         limiter_string = 'not limited';
%     case 2 % normal limiting
%         limiter_string = 'limited';
%     otherwise
%         error('Invalid limiting option');
% end
% FCT_legend_string = ['FCT, ',limiter_string];
% legend_entries = char(legend_entries,FCT_legend_string);
% legend(legend_entries);