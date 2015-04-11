close all; clear; clc;

%% User Options
%--------------------------------------------------------------------------
% mesh options
%--------------------------------------------------------------------------
nel = 10; % number of elements
%--------------------------------------------------------------------------
% method options
%--------------------------------------------------------------------------
% low_order_scheme: 1 = algebraic low-order scheme
%                   2 = graph-theoretic low-order scheme
% high_order_scheme: 1 = Galerkin
%                    2 = Entropy viscosity
%
low_order_scheme  = 2; % low-order scheme option
high_order_scheme = 2; % high-order scheme option
cE = 1.0; % coefficient for entropy residual in entropy viscosity
cJ = 1.0; % coefficient for jumps in entropy viscosity
entropy       = @(u) 0.5*u.^2; % entropy function
entropy_deriv = @(u) u;        % derivative of entropy function
impose_DirichletBC_strongly = true; % impose Dirichlet BC strongly?
%--------------------------------------------------------------------------
% quadrature options
%--------------------------------------------------------------------------
nq = 3; % number of quadrature points
%--------------------------------------------------------------------------
% time options
%--------------------------------------------------------------------------
% temporal_scheme: 1 = explicit Euler
%                  2 = SSPRK(3,3) (Shu-Osher)
%
temporal_scheme = 1; % temporal discretization scheme
CFL = 0.8;       % CFL number
ss_tol = 1.0e-5; % steady-state tolerance
t_end = 1.0;    % max time to run
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP_option: 1 = use DMP
%             2 = max/min(DMP,CMP)
% limiting_option: 0 = set limiting coefficients to 0 (no correction)
%                  1 = set limiting coefficients to 1 (full correction)
%                  2 = use Zalesak's limiter
%                  3 = use Josh's limiter
%
DMP_option = 2;       % DMP option
limiting_option = 2;  % limiter option
%--------------------------------------------------------------------------
% plot options
%--------------------------------------------------------------------------
compute_low_order  = true; % compute and plot low-order solution?
compute_high_order = true; % compute and plot high-order solution?
compute_FCT        = false; % compute and plot FCT solution?
plot_viscosity     = false; % plot viscosities?
plot_low_order_transient  = false; % plot low-order transient?
plot_high_order_transient = true; % plot high-order transient?
plot_FCT_transient        = false; % plot FCT transient?
pausetime = 0.1;                   % time to pause for transient plots
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: pure absorber without source
%            2: void without source -> absorber without source
% IC_option: 0: zero
%            1: exponential pulse
%            2: exponential and square pulse
%
problemID = 3; % problem ID (if any)
periodic_BC = false; % option for periodic BC; otherwise Dirichlet
IC_option = 0; % initial solution option
x_min = 0.0;   % left end of domain
x_max = 1.0;   % right end of domain
mu     = 1;    % cos(angle)
sigma  = @(x) 1;   % cross section function
source = @(x,t) 0; % source function
inc    = 0;    % incoming flux
speed  = 1;    % advection speed
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
save_exact_solution      = false; % option to save exact solution 
save_low_order_solution  = false; % option to save low-order solution
save_high_order_solution = false; % option to save high-order solution
save_FCT_solution        = false; % option to save FCT solution
%-------------------------------------------------------------------------

%% Setup

% set parameters if a problem ID was chosen
switch problemID
    case 0 % custom
        % do nothing, parameters in input section are used
    case 1 % pure absorber without source
        x_min = 0.0;
        x_max = 1.0;
        periodic_BC = false;
        inc    = 1.0;
        mu     = 1.0;
        sigma  = @(x) 10.0;
        source = @(x,t) 0.0;
        speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        exact = @(x,t) (t>=x).*exp(-10*x);
    case 2 % void without source -> absorber without source
        x_min = 0.0;
        x_max = 1.0;
        periodic_BC = false;
        inc    = 1.0;
        mu     = 1.0;
        sigma  = @(x) 10.0*(x >= 0.5);
        source = @(x,t) 0.0;
        speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        exact = @(x,t) (t>=x).*((x<0.5) + (x>=0.5).*(exp(-10*(x-0.5))));
    case 3 % void with source -> absorber without source
        x_min = 0.0;
        x_max = 1.0;
        periodic_BC = false;
        inc    = 0.0;
        mu     = 1.0;
        sigma  = @(x) 10.0*(x >= 0.5);
        source = @(x,t) 1.0*(x < 0.5);
        speed  = 1;
        exact_solution_known = true;
        exact = @(x,t) x.*(x<0.5) + 0.5*exp(-10*(x-0.5)).*(x>=0.5);
    case 5 % MMS-1
        x_min = 0.0;
        x_max = 1.0;
        periodic_BC = false;
        inc    = 0.0;
        mu     = 1.0;
        sigma  = @(x) 1.0;
        source = @(x,t) sin(pi*x)+pi*t*cos(pi*x)+t*sin(pi*x);
        speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        exact = @(x,t) t*sin(pi*x);
    otherwise
        error('Invalid problem ID chosen');
end

% compute mesh quantities
% number of dofs
if periodic_BC
    n_dof = nel;
else
    n_dof = nel + 1;
end

% mesh
x = linspace(x_min,x_max,nel+1)'; % positions of dofs/cell vertices
dx = diff(x);                     % element sizes
if periodic_BC
    x(end)=[];
end

% set initial conditions function and compute initial solution
switch IC_option
    case 0 % zero
        IC = @(x) zeros(length(x),1);
    case 1 % exponential pulse
        IC = @(x) exp(-200*(x-0.5).^2);
    case 2 % exponential and square pulse
        IC = @(x) exp(-200*(x-0.3).^2)+(x>=0.6).*(x<=0.8);
    case 3 % square pulse
        IC = @(x) (x>=0.4).*(x<0.7);
    otherwise
        error('Invalid IC option');
end
% compute initial solution
u0 = IC(x);

% get quadrature points and weights and evaluate basis functions
[zq,wq]  = get_GL_quadrature(nq);
[v,dvdz] = get_lagrange_basis(zq,2); 
Jac = 0.5*dx; % Jacobians of reference cell transformations

%% Material properties

% compute min and max sigma and source in the support of i
sigma_min  = 1e15*ones(n_dof,1);
sigma_max  = zeros(n_dof,1);
source_min = 1e15*ones(n_dof,1);
source_max = zeros(n_dof,1);
for iel = 1:nel
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    
    % compute cross section and source at each quadrature point
    sigma_cell = sigma(xq);
    source_cell = source(xq,0);
    
    % compute max and min on cell
    sigma_cell_min = min(sigma_cell);
    sigma_cell_max = max(sigma_cell);
    source_cell_min = min(source_cell);
    source_cell_max = max(source_cell);
    
    % update max/min for each dof on cell
    sigma_min(iel:iel+1) = min(sigma_min(iel:iel+1),sigma_cell_min);
    sigma_max(iel:iel+1) = min(sigma_max(iel:iel+1),sigma_cell_max);
    source_min(iel:iel+1) = min(source_min(iel:iel+1),source_cell_min);
    source_max(iel:iel+1) = min(source_max(iel:iel+1),source_cell_max);
end

%% Assembly

% create connectivity array
connectivity = [linspace(1,nel,nel)' linspace(2,nel+1,nel)'];
if periodic_BC
    connectivity(end,:) = [connectivity(end,1) 1];
end

% determine if linear system needs to be modified for Dirichlet BC
modify_for_weak_DirichletBC = ~periodic_BC && ~impose_DirichletBC_strongly;
modify_for_strong_DirichletBC = ~periodic_BC && impose_DirichletBC_strongly;

% assemble mass matrices, inviscid and low-order steady-state matrices,
% low-order viscosity and corresponding artificial diffusion matrix
[MC,ML,A,AL,DL,viscL] = build_matrices(...
    nq,zq,wq,v,dvdz,Jac,x,dx,nel,n_dof,connectivity,...
    mu,speed,sigma,low_order_scheme,modify_for_weak_DirichletBC);

% assemble steady-state rhs at time 0
b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
    source,mu,speed,0,inc,modify_for_weak_DirichletBC);

% check that speed is nonzero; otherwise, infinite dt will be computed
if (speed == 0)
    error('Speed cannot be zero; infinite dt will be computed');
end

% compute dt using CFL condition, even for implicit, just so that time
% steps will be equal between explicit and implicit, for comparison
max_speed_dx = 0.0;
for i = 1:n_dof
    max_speed_dx = max(max_speed_dx, AL(i,i)/ML(i,i));
end
dtCFL = CFL / max_speed_dx;

% compute transient sytem matrices and modify for Dirichlet BC
ALtr = ML;
ALtr_mod = ALtr;
if modify_for_strong_DirichletBC
    ALtr_mod(1,:)=0; ALtr_mod(1,1)=1;
end

%% Low-order Solution

if (compute_low_order)
    fprintf('\nComputing low-order solution...\n\n');
    u_old = u0;
    t = 0;
    time_step = 0;
    reached_end_of_transient = false;
    while (~reached_end_of_transient)
        % advance time step index
        time_step = time_step+1;
        
        % shorten time step if necessary
        if (t+dtCFL >= t_end)
            dt = t_end - t;
            reached_end_of_transient = true;
        else
            dt = dtCFL;
            ALtr = ML;
            ALtr_mod = ALtr;
            if modify_for_strong_DirichletBC
                ALtr_mod(1,:)=0; ALtr_mod(1,1)=1;
            end
            reached_end_of_transient = false;
        end
        fprintf('Time step %i: t = %f->%f',time_step,t,t+dt);
        
        % perform step
        switch temporal_scheme
            case 1 % explicit Euler
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                uL = low_order_step(u_old,AL,ALtr_mod,ML,b,dt,0,inc,...
                    modify_for_strong_DirichletBC);
            case 2 % SSP3
                % stage 1
                u_old_stage = u_old;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                uL_stage = low_order_step(u_old_stage,AL,ALtr_mod,ML,b,...
                    dt,0,inc,modify_for_strong_DirichletBC);
                % stage 2
                u_old_stage = uL_stage;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+dt,inc,modify_for_weak_DirichletBC);
                uL_stage = low_order_step(u_old_stage,AL,ALtr_mod,ML,b,...
                    dt,0,inc,modify_for_strong_DirichletBC);
                % stage 3
                u_old_stage = 0.75*u_old + 0.25*uL_stage;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+0.5*dt,inc,modify_for_weak_DirichletBC);
                uL_stage = low_order_step(u_old_stage,AL,ALtr_mod,ML,b,...
                    dt,0,inc,modify_for_strong_DirichletBC);
                % final combination
                uL = 1/3*u_old + 2/3*uL_stage;
            otherwise
                error('Invalid temporal discretization scheme');
        end
        
        % test steady-state convergence
        ss_err = norm(uL-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        reached_steady_state = ss_err < ss_tol;
        if (reached_steady_state)
            break;
        end
        
        % plot
        if (plot_low_order_transient)
            plot(x,uL);
            legend('Low-order','Location','Best');
            axis([0 len 0 1]);
            pause(pausetime);
        end
        
        % reset u_old
        u_old = uL;
        t = t+dt;
    end
    % save uL for plotting because "uL" is used in FCT loop
    uL_final = uL;
end

%% High-order solution

if (compute_high_order)
    fprintf('\nComputing high-order solution...\n\n');
    u_old   = u0;
    u_older = u0;
    t = 0;
    time_step = 0;
    dt_old = dtCFL;
    reached_end_of_transient = false;
    while (~reached_end_of_transient)
        % advance time step index
        time_step = time_step+1;
        
        % shorten time step if necessary
        if (t+dtCFL >= t_end)
            dt = t_end - t;
            reached_end_of_transient = true;
        else
            dt = dtCFL;
            reached_end_of_transient = false;
        end
        fprintf('Time step %i: t = %f->%f',time_step,t,t+dt);
        
        switch temporal_scheme
            case 1 % explicit Euler
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                [uH,DH] = high_order_step(u_older,u_old,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
            case 2 % SSP3
                % stage 1
                u_old_stage = u_old;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                % stage 2
                u_old_stage = uH_stage;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+dt,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                % stage 3
                u_old_stage = 0.75*u_old + 0.25*uH_stage;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+0.5*dt,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                % final combination
                uH = 1/3*u_old + 2/3*uH_stage;
            otherwise
                error('Invalid temporal discretization scheme');
        end
        
        % test steady-state convergence
        ss_err = norm(uH-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        reached_steady_state = ss_err < ss_tol;
        if (reached_steady_state)
            break;
        end
        
        % plot
        if (plot_high_order_transient)
            plot(x,uH);
            axis([-inf inf -inf 0.5]);
            legend('High-order');
            pause(pausetime);
        end
        
        % reset u_old and u_older and advance time
        u_older = u_old;
        u_old   = uH;
        dt_old = dt;
        t = t+dt;
    end
    
    % save uH for plotting because "uH" is used in FCT loop
    uH_final = uH;
end

%% FCT solution

if (compute_FCT)
    fprintf('\nComputing FCT solution...\n\n');
    u_old   = u0;
    u_older = u0;
    t = 0;
    time_step = 0;
    dt_old = dtCFL;
    reached_end_of_transient = false;
    while (~reached_end_of_transient)
        % advance time step index
        time_step = time_step+1;
        
        % shorten time step if necessary
        if (t+dtCFL >= t_end)
            dt = t_end - t;
            reached_end_of_transient = true;
        else
            dt = dtCFL;
            reached_end_of_transient = false;
        end
        fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
        
        switch temporal_scheme
            case 1 % explicit Euler
                % perform high-order step
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                [uH,DH] = high_order_step(u_older,u_old,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);

                % perform FCT step
                uFCT = FCT_step_explicit(u_old,uH,dt,...
                    ML,MC,AL,DH,DL,b,inc,speed,sigma_min,sigma_max,...
                    source_min,source_max,DMP_option,limiting_option,...
                    periodic_BC,impose_DirichletBC_strongly);
            case 2 % SSP3
                % stage 1
                u_old_stage = u_old;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                    ML,MC,AL,DH,DL,b,inc,speed,sigma_min,sigma_max,...
                    source_min,source_max,DMP_option,limiting_option,...
                    periodic_BC,impose_DirichletBC_strongly);
                
                % stage 2
                u_old_stage = uFCT_stage;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+dt,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                    ML,MC,AL,DH,DL,b,inc,speed,sigma_min,sigma_max,...
                    source_min,source_max,DMP_option,limiting_option,...
                    periodic_BC,impose_DirichletBC_strongly);
                
                % stage 3
                u_old_stage = 0.75*u_old + 0.25*uFCT_stage;
                u_older_stage = u_older;
                b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
                    source,mu,speed,t+0.5*dt,inc,modify_for_weak_DirichletBC);
                [uH_stage,DH] = high_order_step(u_older_stage,u_old_stage,viscL,dx,...
                    x,mu,sigma,source,inc,dt,dt_old,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
                    A,b,MC,0,time_step,high_order_scheme,periodic_BC,impose_DirichletBC_strongly);
                uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                    ML,MC,AL,DH,DL,b,inc,speed,sigma_min,sigma_max,...
                    source_min,source_max,DMP_option,limiting_option,...
                    periodic_BC,impose_DirichletBC_strongly);
                
                % final combination
                uFCT = 1/3*u_old + 2/3*uFCT_stage;
            otherwise
                error('Invalid temporal discretization scheme');
        end
        
        % test steady-state convergence
        ss_err = norm(uFCT-u_old,2);
        fprintf('\tnorm(u_new - u_old) = %e\n',ss_err);
        reached_steady_state = ss_err < ss_tol;
        if (reached_steady_state)
            break;
        end
        
        % plot
        if (plot_FCT_transient)
            plot(x,uFCT,'g');
            legend('FCT','Location','Best');
            pause(pausetime);
        end
        
        % reset u_old and u_older and advance time
        u_older = u_old;
        u_old   = uFCT;
        dt_old = dt;
        t = t+dt;
    end
end

%% Plot Solution

figure(1); clf; hold on;

% plot exact solution
if (exact_solution_known)
    xx = linspace(x_min,x_max,1000)';
    u_exact = exact(xx,t_end);
    plot(xx,u_exact,'k-');
    legend_entries = char('Exact');
end

% plot low-order solution
if (compute_low_order)
    plot(x,uL_final,'r-s');
    if (exact_solution_known)
        legend_entries = char(legend_entries,'Low-order');
    else
        legend_entries = char('Low-order');
    end
end

% plot high-order solution
if (compute_high_order)
    plot(x,uH_final,'b-+');
    legend_entries = char(legend_entries,'High-order');
end

% plot FCT solution
if (compute_FCT)
    plot(x,uFCT,'g-x');
    switch limiting_option
        case 0 % no correction
            limiter_string = 'no correction';
        case 1 % full correction (no limiting)
            limiter_string = 'not limited';
        case 2 % Zalesak limiter
            limiter_string = 'Zalesak limiter';
        case 3 % Josh limiter
            limiter_string = 'Josh limiter';
        otherwise
            error('Invalid limiting option');
    end
    FCT_legend_string = ['FCT, ',limiter_string];
    legend_entries = char(legend_entries,FCT_legend_string);
end

% legend
legend(legend_entries,'Location','Best');

%% Plot viscosity

% plot viscosity if requested
if (plot_viscosity)
    figure(2); clf;
    
    % plot low-order viscosity if available
    if low_order_scheme == 2 || high_order_scheme == 2
        semilogy(x_center,viscL);
        legend_entries = char('Low-order viscosity');
    end
    
    hold on;
    
    % plot high-order viscosity if available
    if (high_order_scheme == 2)
        semilogy(x_center,viscE,'x');
        legend_entries = char(legend_entries,'Entropy viscosity');
    end
    
    % legend
    legend(legend_entries,'Location','Best');
end

%% Output

% save exact solution
if (save_exact_solution)
    exact_file = 'output/uexact.txt';
    dlmwrite(exact_file,[xx,u_exact],' ');
end

% determine string to be appended to results for time discretization
switch temporal_scheme
    case 1
        time_string = 'FE';
    case 2
        time_string = 'SSPRK33';
    otherwise
        error('Invalid time discretization scheme');
end

% save low-order solution
if (save_low_order_solution)
    if (compute_low_order)
        low_order_file = ['output/uL_',time_string,'.txt'];
        dlmwrite(low_order_file,[x,uL_final],' ');
    end
end

% determine string to be appended for high-order scheme
switch high_order_scheme
    case 1
        high_order_string = 'Gal';
    case 2
        high_order_string = 'EV';
    otherwise
        error('Invalid high-order scheme chosen');
end

% save high-order solution
if (save_high_order_solution)
    if (compute_high_order)
        high_order_file = ['output/uH_',high_order_string,'_',time_string,'.txt'];
        dlmwrite(high_order_file,[x,uH_final],' ');
    end
end

% save FCT solution
if (save_FCT_solution)
    if (compute_FCT)
        FCT_file = ['output/uFCT_',high_order_string,'_',time_string,'.txt'];
        dlmwrite(FCT_file,[x,uFCT],' ');
    end
end