close all; clc; clear;
%% Introduction

% Solves a transient transport equation:
%   MC*du/dt-(K  )*u = b (high-order system)
%   ML*du/dt-(K+D)*u = b (low-order system)
% using the theta-scheme for time discretization

%% User Options

nel_init = 50; % number of elements in first cycle
n_cycle  = 1;  % number of refinement cycles

problemID = 1; % 1 = exponential and square wave
               
theta = 0.0; % 0:   explicit euler
             % 1:   implicit euler
             % 1/2: crank nicholson
CFL = 0.8;   % CFL number
t_end = 0.5; % end time
force_dt = false; % option to use dt_user instead of computing from CFL
dt_user = 0.001;  % dt used if force_dt == true

limiting_option = 2; % 0 = set limiting coefficients to 0 (no correction)
                     % 1 = set limiting coefficients to 1 (full correction)
                     % 2 = compute limiting coefficients normally

save_convergence_data = 0; % option to save convergence data (0=no,1=yes)

%% Setup

% determine problem parameters from ID
switch problemID
    case 1
        len=1;          % length of domain
        omega=1;        % omega
        sigma = @(x) 0; % sigma
        src=0;          % source
        inc=1;          % incoming flux
        speed=1;        % transport speed
        IC = @IC_exponential_and_square;
        %IC = @IC_square;
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
    n_dof = nel;                  % number of dofs
    h(cycle) = len/nel;           % mesh size
    fprintf('Cycle %i: n_dof = %i\n',cycle,n_dof);
    
    % build matrices
    [MC,ML,K,D,b] = build_matrices(len,nel,omega,sigma,src,speed);
    ML = diag(ML);
    AL = -(K+D);    % low-order steady-state system matrix
    AH = -K;        % high-order steady-state system matrix
    
    % compute dt using CFL condition for explicit Euler, even if not
    % using explicit Euler, just for consistency in dt for comparison
    dt0 = abs(CFL*ML(1,1)/AL(1,1));
    for i = 2:n_dof
        dt0 = min(dt0, abs(CFL*ML(i,i)/AL(i,i)));
    end
    % if user provided dt, use that dt instead of that computed from CFL
    if (force_dt)
        dt0 = dt_user;
    end
    
    % initial conditions for pseudotransient (equal to exact steady-state solution)
    x = linspace(0,len-h(cycle),n_dof);
    u0 = IC(x);
    
    % low-order solve
    %======================================================================
    fprintf('\tComputing low-order solution...\n');
    u_old = u0;
    t = 0;
    in_transient = true;
    while (in_transient)
        if (t+dt0 >= t_end)
            dt = t_end - t;
            in_transient = false;
        else
            dt = dt0;
        end
        fprintf('\t\tTime step t = %f->%f\n',t,t+dt);
        % compute system matrix
        ALtr = ML + theta*dt*AL;
        % compute rhs
        rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b;
        % solve modified system
        uL = ALtr \ rhs;
        % reset u_old
        u_old = uL;
        t = t+dt;
    end

    % high-order solve
    %======================================================================
    fprintf('\tComputing high-order solution...\n');
    u_old = u0;
    t = 0;
    in_transient = true;
    while (in_transient)
        if (t+dt0 >= t_end)
            dt = t_end - t;
            in_transient = false;
        else
            dt = dt0;
        end
        fprintf('\t\tTime step t = %f->%f\n',t,t+dt);
        % compute system matrix
        AHtr = MC + theta*dt*AH;
        % compute rhs
        rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
        % solve system
        uH = AHtr \ rhs;
        % reset u_old
        u_old = uH;
        t = t+dt;
    end
    
    % FCT solve
    %======================================================================
    fprintf('\tComputing FCT solution...\n');
    u_old = u0;
    uFCT  = u0;
    t = 0;
    in_transient = true;
    while (in_transient)
        if (t+dt0 >= t_end)
            dt = t_end - t;
            in_transient = false;
        else
            dt = dt0;
        end
        fprintf('\t\tTime step t = %f->%f\n',t,t+dt);
        
        % compute system matrices
        ALtr = ML + theta*dt*AL;
        AHtr = MC + theta*dt*AH;
        
        % compute auxiliary solution
        %----------------
        % compute rhs
        rhs = (ML - (1-theta)*dt*AL)*u_old + (1-theta)*dt*b;
        % solve for aux. solution
        u_aux = ML \ rhs;
    
        % high-order solve
        %----------------
        % compute rhs
        rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
        % solve modified system
        uH_FCT_loop = AHtr \ rhs;
        
        % FCT solve
        %----------------
        % compute flux correction matrix
        F = flux_correction_matrix(u_old,uH_FCT_loop,dt,D,MC,theta);
        % cancel down gradient as Kuzmin suggested
        %F = cancel_down_gradient(F,u_aux);
        % compute limiting coefficients
        switch limiting_option
            case 0 % no correction
                Flim = 0*F;
            case 1 % full correction (no limiting)
                Flim = F;
            case 2 % normal limiting
                Flim = limit_fluxes(F,u_aux,ML);
            otherwise
                error('Invalid limiting option');
        end
        % enforce LED as Kuzmin suggested
        %Flim = enforce_LED(Flim,u_aux);
        % compute correction rhs
        rhs = ML*u_aux + theta*dt*b + sum(Flim,2);
        % solve for FCT solution
        uFCT = ALtr \ rhs;

        % plot
        figure(1); clf; hold on;
        plot(x,uFCT);
        
        % reset u_old
        u_old = uFCT;
        t = t+dt;
    end

%     % compute error
%     %======================================================================
%     uL_err(cycle)   = compute_error(uL,x);
%     uH_err(cycle)   = compute_error(uH,x);
%     uFCT_err(cycle) = compute_error(uFCT,x);
end

%% Convergence Rates
% 
% fprintf('\nLow-order Convergence:\n')
% for cycle = 1:n_cycle
%     if (cycle == 1)
%         fprintf('Cycle %i: L2 error = %e\n',cycle,uL_err(cycle));
%     else
%         rate = log(uL_err(cycle)/uL_err(cycle-1))/log(h(cycle)/h(cycle-1));
%         fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uL_err(cycle),rate);
%     end
% end
% fprintf('\nHigh-order Convergence:\n')
% for cycle = 1:n_cycle
%     if (cycle == 1)
%         fprintf('Cycle %i: L2 error = %e\n',cycle,uH_err(cycle));
%     else
%         rate = log(uH_err(cycle)/uH_err(cycle-1))/log(h(cycle)/h(cycle-1));
%         fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uH_err(cycle),rate);
%     end
% end
% fprintf('\nFCT Convergence:\n')
% for cycle = 1:n_cycle
%     if (cycle == 1)
%         fprintf('Cycle %i: L2 error = %e\n',cycle,uFCT_err(cycle));
%     else
%         rate = log(uFCT_err(cycle)/uFCT_err(cycle-1))/log(h(cycle)/h(cycle-1));
%         fprintf('Cycle %i: L2 error = %e, rate = %f\n',cycle,uFCT_err(cycle),rate);
%     end
% end
% 
% if (save_convergence_data)
%     % save convergence data to file
%     dlmwrite('output/convergence.csv',[h,uL_err,uH_err,uFCT_err]);
% end

%% Plot

figure(1); clf; hold on;
% exact solution
x_exact = linspace(0,len,1000);
u_exact = exact_solution(IC,x_exact,t_end,len);
plot(x_exact,u_exact,'k-');
% numerical solutions
plot(x,uL,'r-s');
% plot(x,uH,'b-+');
plot(x,uFCT,'g-x');
% plot legend
legend_entries = char('Exact');
legend_entries = char(legend_entries,'Low-order');
% legend_entries = char(legend_entries,'High-order');
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
legend(legend_entries,'Location','Best');